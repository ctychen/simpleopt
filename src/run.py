import time
import sys
import numpy as np
import os

import pandas as pd
import scipy

from multiprocessing import Process

try:
    runMode = os.environ["runMode"]
except:
    runMode = 'local'
    os.environ["runMode"] = runMode

if (runMode == 'docker'):
    FreeCADPath = '/usr/lib/freecad-daily/lib'
    BlenderPath = '/usr/lib/blender'
    OpenSCADPath = '/usr/bin/openscad'
else:
    FreeCADPath = '/Applications/FreeCAD.app/Contents/Resources/lib'
    BlenderPath = '/Applications/Blender.app/Contents/MacOS/blender'

sys.path.append(FreeCADPath)
sys.path.append(BlenderPath)
sys.path.append(OpenSCADPath)
print(sys.path)

import trimesh

import CADClass

import Solid
import ForwardModel
import OptModel
import ObjectiveFunctionTools

import multiprocessing

class RunSetup_MeshHF:
    def __init__(self):
        g_obj = lambda qvals: max(qvals) #+ qvals.count(max(qvals)) #maybe changing obj function helps??

        stpPath = "unit_test_cube.step" #"unit_test_cone.step" 
        # stpPath = "test_sphere.step"
        # stpPath = "test_pfc_block.step" #"unit_test_pfc.step" #for multiple directions 

        stlPath = " " #"box.stl"

        # qDirIn = [[0.707, -0.707, 0.0], [-0.707, -0.707, 0.0]] #[m]
        qDirIn = [0.0, -1.0, 0.0] #[m]
        # qDirIn = [0.707, -0.707, 0.0] #[m]
        # qDirIn = [-0.707, -0.707, 0.0] #[m]
        # qDirIn = [0.0, -0.707, 0.707] #[m]
        # qDirIn = [0.0, -0.707, -0.707] #[m]
        qMagIn = 10.0 #[W/m^2]
        # qMagIn = [10.0, 10.0] #[W/m^2]
#
        self.box = Solid.MeshSolid(stlPath, stpPath)

        self.fwd = ForwardModel.ForwardModel_MeshHF(self.box, qMagIn, qDirIn, hfMode='uniform') 
        # self.fwd = ForwardModel.ForwardModel_MeshHF(self.box, qMagIn, qDirIn, hfMode='uniform_multiple') 
        # self.fwd = ForwardModel.ForwardModel_MeshHF(self.box, qMagIn, qDirIn, hfMode='exponnorm') 
        self.opt = OptModel.OptModel_MeshHF()

        print(f"Number of cores: {multiprocessing.cpu_count()}")

        return
    

    def makeDirectories(self, runID, coefficientsList):
        c0 = coefficientsList[0]
        c1 = coefficientsList[1]
        c2 = coefficientsList[2]
        c3 = coefficientsList[3]
        c4 = coefficientsList[4]
        runName = runID + f'_c0_{c0:.5f}_c1_{c1:.5f}_c2_{c2:.5f}_c3_{c3:.5f}_c4_{c4:.5f}' #f'_c0_{c0:.2f}_c1_{c1:.2f}_c2_{c2:.2f}_c3_{c3:.2f}'
        runName = runName.replace(".", "-")

        directoryName = f"{runName}" #running this within docker container means can't save to external without bindmount aaa

        os.makedirs(directoryName)

        print(f"Made directory: {directoryName}")

        return directoryName

    def makeSweepCSV(self, c0List, c1List, c2List, c3List, maxhfList, runID):
        data = {
            'c0': c0List,
            'c1': c1List,
            'c2': c2List,
            'c3': c3List,
            'maxhf': maxhfList
        }
        df = pd.DataFrame(data)
        df.to_csv(f"{runID}_results.csv", index=None, header=True)

        return 
    

    def runOptimization(self):

        self.box.processSolid()
        trimeshSolid = self.box.trimeshSolid

        #find which faces are on corners - don't want to take those into account for difference in normals? 
        #reasoning is more than 1 adjacent face, and more than 1 normal to compare as a result? 
        #maybe this is not the way to go? 
        # corner_faces = self.findCornerFaces(self.box.trimeshSolid)
        numVertices = len(trimeshSolid.vertices) 
        print(f"Shape of trimesh vertices: {trimeshSolid.vertices.shape}")  
        numFaces = len(trimeshSolid.faces)
        print(f"Number of trimesh solid vertices: {numVertices}")
        print(f"Number of trimesh solid faces: {numFaces}")

        def facesToKeep(trimeshSolid):
            #this should be for faces to NOT use (so they shouldn't move)
            return []


        numFaces = len(trimeshSolid.faces)  #for normalizing terms in objective function

        def sweep_coefficients_and_record_output(coefficients_list, idx_to_vary, sweep_values):
            """
            run through possible values of one coefficient while keeping the rest fixed
            do optimization and save results and attempted combinations to csv
            """

            c0_runvals = []
            c1_runvals = []
            c2_runvals = []
            c3_runvals = []
            c4_runvals = []
            maxhf_vals = []

            for val in sweep_values:
                my_trimeshSolid = trimeshSolid.copy()
                coefficients_list[idx_to_vary] = val
                # directoryName = self.makeDirectories(f"vertexscans1_0mm/imc_and_maxvtx/c2_c3_test_", coefficients_list)
                directoryName = self.makeDirectories(f"1_0mm_sphere_testing/imc_and_vtxdefects_vtxnotnormalized/test_", coefficients_list)
                #meshHFOpt(self, hfObjectiveFcn, constraint, updateHFProfile, calcHFAllMesh, calcMaxHF, calcEnergy, meshObj, coefficientsList, threshold, delta, id):
                self.opt.meshHFOpt(
                    facesToKeep,
                    self.fwd.makeHFProfile,
                    my_trimeshSolid, 
                    coefficients_list,
                    threshold=0.00000001, 
                    delta=0.01, 
                    fwdModel=self.fwd,
                    id=directoryName
                )
                c0_runvals.append(coefficients_list[0])
                c1_runvals.append(coefficients_list[1])
                c2_runvals.append(coefficients_list[2])
                c3_runvals.append(coefficients_list[3])
                c4_runvals.append(coefficients_list[4])
                # maxhf_vals.append(maxHF)
            
            # self.makeSweepCSV(c0_runvals, c1_runvals, c2_runvals, c3_runvals, c4_runvals, f"newspheretest{idx_to_vary}/1/{idx_to_vary}")
            return 
        
        def big_sweep_coefficients_and_record_output(sweep_values_c0, sweep_values_c1, sweep_values_c2, sweep_values_c3, sweep_values_c4, id):
            """
            run through possible values of one coefficient while keeping the rest fixed
            do optimization and save results and attempted combinations to csv
            """

            c0_runvals = []
            c1_runvals = []
            c2_runvals = []
            c3_runvals = []
            c4_runvals = []
            maxhf_vals = []

            # for c0 in sweep_values_c0:
                # for c1 in sweep_values_c1: 
                #     for c2 in sweep_values_c2:
                #         for c3 in sweep_values_c3: 
                #             my_trimeshSolid = trimeshSolid.copy()
                #             # coefficients_list[0] = c0
                #             coefficients_list = [c0, c1, c2, c3, 0.0]
                #             directoryName = self.makeDirectories(f"randomspheretest_{id}/test_", coefficients_list)
                            
                #             self.opt.meshHFOpt(
                #                 objectiveFunction,  
                #                 facesToKeep,
                #                 self.fwd.makeHFProfile,
                #                 self.fwd.calculateAllHF,
                #                 self.fwd.filteredCalculateMaxHF, #self.fwd.calculateMaxHF,
                #                 self.fwd.calculateIntegratedEnergy,
                #                 my_trimeshSolid, 
                #                 coefficients_list,
                #                 threshold=0.00000001, 
                #                 delta=0.01, 
                #                 id=directoryName
                #             )
                #             c0_runvals.append(coefficients_list[0])
                #             c1_runvals.append(coefficients_list[1])
                #             c2_runvals.append(coefficients_list[2])
                #             c3_runvals.append(coefficients_list[3])
                            # maxhf_vals.append(maxHF)

            for c2 in sweep_values_c2:
                for c4 in sweep_values_c4: 
                    my_trimeshSolid = trimeshSolid.copy()
                    # coefficients_list[0] = c0
                    coefficients_list = [0, 0, c2, 0, c4]
                    directoryName = self.makeDirectories(f"randomspheretest2_{id}/test_", coefficients_list)
                    
                    self.opt.meshHFOpt(
                        facesToKeep,
                        self.fwd.makeHFProfile,
                        my_trimeshSolid, 
                        coefficients_list,
                        threshold=0.00000001, 
                        delta=0.01, 
                        fwdModel=self.fwd,
                        id=directoryName
                    )
                    c0_runvals.append(coefficients_list[0])
                    c1_runvals.append(coefficients_list[1])
                    c2_runvals.append(coefficients_list[2])
                    c3_runvals.append(coefficients_list[3])
                    c4_runvals.append(coefficients_list[4])
                    # maxhf_vals.append(maxHF)
            
            self.makeSweepCSV(c0_runvals, c1_runvals, c2_runvals, c3_runvals, c4_runvals, maxhf_vals, f"randomspheretest2_{id}/{id}")
            return 


        coefficients_list = [0, 0, 6.67 * 2 *  94.24777960769379, 10, 0] #[0, 0, 0, 20, 10]
        # import random
        # sweep_c2 = random.uniform(1000, 10000) 
        sweep_c2 = [6.67 * 2 *  94.24777960769379]
        sweep_coefficients_and_record_output(coefficients_list, 2, sweep_c2)    
        return 


if __name__ == '__main__':

    t0 = time.time()

    setup = RunSetup_MeshHF()

    """
    Run cube mesh operation test
    """
    setup.runOptimization()

    print(f"Time elapsed: {time.time() - t0}")





        
