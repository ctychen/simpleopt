import time
import sys
import numpy as np
import os

import pandas as pd

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import plotly.express as px

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

    # HEATPath = '/root/source/HEAT'
else:
    #FreeCADPath = '/Applications/FreeCAD.app' #'/usr/lib/freecad-python3/lib'
    FreeCADPath = '/Applications/FreeCAD.app/Contents/Resources/lib'
    BlenderPath = '/Applications/Blender.app/Contents/MacOS/blender'
    # HEATPath = '/Users/cchen/Desktop/HEAT'

sys.path.append(FreeCADPath)
sys.path.append(BlenderPath)
sys.path.append(OpenSCADPath)
print(sys.path)

import trimesh

import CADClass

import Solid
import ForwardModel
import OptModel

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
        #self.box = Solid.MeshSolid(stlPath, stpPath) #normally, use this one!

        self.fwd = ForwardModel.ForwardModel_MeshHF(self.box, qMagIn, qDirIn, hfMode='uniform') 
        # self.fwd = ForwardModel.ForwardModel_MeshHF(self.box, qMagIn, qDirIn, hfMode='uniform_multiple') 
        # self.fwd = ForwardModel.ForwardModel_MeshHF(self.box, qMagIn, qDirIn, hfMode='exponnorm') 
        self.opt = OptModel.OptModel_MeshHF()

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
        # os.makedirs(f"{directoryName}/images")

        print(f"Made directory: {directoryName}")

        #where we are saving generated vtk's
        # directoryName = f"{directoryName}/vtks"
        # os.makedirs(directoryName)

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
        adjacency = trimeshSolid.face_adjacency

        def calculateNormalsDiff(trimeshSolid):
            """
            calculate differences between normal vectors of adjacent faces - ideally this should also be minimized?  
            also find max diff btwn adjacent normals
            """
            #"(n, 2) list of face indices. Each pair of faces in the list shares an edge, making them adjacent"
            # adjacency = trimeshSolid.face_adjacency 

            normals = trimeshSolid.face_normals
            # normalsDiff = normals[adjacency[:, 0]] - normals[adjacency[:, 1]]
            # normalsDiffMagnitude = np.linalg.norm(normalsDiff, axis=1)

            vertex_defects = trimeshSolid.vertex_defects
            sumVertexDefects = np.sum(np.abs(vertex_defects))
            maxVertexDefect = np.max(np.abs(vertex_defects))    

            normals_0 = normals[adjacency[:, 0]]
            normals_1 = normals[adjacency[:, 1]]
            dot_product = np.einsum('ij,ij->i', normals_0, normals_1)
            clipped_dot_product = np.clip(dot_product, -1.0, 1.0)
            allAnglesBetweenNormals = np.arccos(clipped_dot_product)
            maxAngleBetweenNormals = np.max(allAnglesBetweenNormals)

            return sumVertexDefects, maxVertexDefect, maxAngleBetweenNormals 
            #return normalsDiffMagnitude, maxAngleBetweenNormals

        # def calculateNormalsDiff(trimeshSolid):
        #     """
        #     calculate differences between normal vectors of adjacent faces - ideally this should also be minimized?  
        #     also find max diff btwn adjacent normals
        #     """
        #     #"(n, 2) list of face indices. Each pair of faces in the list shares an edge, making them adjacent"
        #     # adjacency = trimeshSolid.face_adjacency 

        #     # normals = trimeshSolid.face_normals
        #     # normalsDiff = normals[adjacency[:, 0]] - normals[adjacency[:, 1]]
        #     # normalsDiffMagnitude = np.linalg.norm(normalsDiff, axis=1)

        #     vertex_defects = trimeshSolid.vertex_defects
        #     sumVertexDefects = np.sum(np.abs(vertex_defects))
        #     maxVertexDefect = np.max(np.abs(vertex_defects))    

        #     # normals_0 = normals[adjacency[:, 0]]
        #     # normals_1 = normals[adjacency[:, 1]]
        #     # dot_product = np.einsum('ij,ij->i', normals_0, normals_1)
        #     # clipped_dot_product = np.clip(dot_product, -1.0, 1.0)
        #     # allAnglesBetweenNormals = np.arccos(clipped_dot_product)
        #     # maxAngleBetweenNormals = np.max(allAnglesBetweenNormals)
        #     maxAngleBetweenNormals = 0

        #     return sumVertexDefects, maxVertexDefect, maxAngleBetweenNormals 
        #     #return normalsDiffMagnitude, maxAngleBetweenNormals
            

        def facesToKeep(trimeshSolid):
            #this should be for faces to NOT use (so they shouldn't move)
            return []

        
        q_mesh_initial = self.fwd.calculateAllHF(trimeshSolid)
        maxHF_initial = self.fwd.filteredCalculateMaxHF(q_mesh_initial, unconstrainedFaces = [])
        sumHF_initial = self.fwd.calculateHFMeshSum(q_mesh_initial)
        # normalsDiff_initial, normalRefDotProducts_initial = calculateNormalsDiff(trimeshSolid)
        # normalsDiff_initial, maxNormalsDiff_initial = calculateNormalsDiff(trimeshSolid)
        sumVertexDefects_initial, maxVertexDefect_initial, maxAngleBetweenNormals_initial = calculateNormalsDiff(trimeshSolid)   
        # normalsPenalty_initial = np.sum(normalsDiff_initial)
        energy_initial = self.fwd.calculateIntegratedEnergy(q_mesh_initial, trimeshSolid)

        numFaces = len(trimeshSolid.faces)  #for normalizing terms in objective function

        print(f"Initial max HF: {maxHF_initial}")
        print(f"Initial sum HF: {sumHF_initial}")
        # print(f"Initial normals penalty: {normalsPenalty_initial}")
        # print(f"Initial max angle between normals: {maxNormalsDiff_initial}")
        print(f"Initial sum vertex defects: {sumVertexDefects_initial}")
        print(f"Initial max vertex defect: {maxVertexDefect_initial}")
        print(f"Initial max angle between normals: {maxAngleBetweenNormals_initial}")
        print(f"Initial energy: {energy_initial}")

        # input()

        def objectiveFunction(trimeshSolid, coefficientsList, unconstrainedFaces):

            t0 = time.time()

            c0 = coefficientsList[0] #max heat flux term
            c1 = coefficientsList[1] #sum heat flux, unweighted, term
            c2 = coefficientsList[2] #normals diff term
            c3 = coefficientsList[3] #max normals term now 

            c4 = coefficientsList[4] #heat flux diff term

            q_mesh_all = self.fwd.calculateAllHF(trimeshSolid)

            # print(f"Time elapsed for q_mesh_all: {time.time() - t0}")

            unconstrainedFaces = [] #only keep this for cases with hf only on top face? 

            maxHFTerm = 0 #c0 * self.fwd.filteredCalculateMaxHF(q_mesh_all, unconstrainedFaces)    #try not dividing by initial value
            sumHFTerm = 0 #c1 * (self.fwd.calculateHFMeshSum(q_mesh_all) / numFaces) 

            # normalsDiff, maxNormalsDiff = calculateNormalsDiff(trimeshSolid)
            # normalsDiffSum = np.sum(normalsDiff)
            # normalsPenalty = c2 * ((normalsDiffSum / normalsPenalty_initial) / numFaces)
            # normalsPenalty = (c2) * (normalsDiffSum / numFaces)

            # sumVertexDefects, maxAngleBetweenNormals = calculateNormalsDiff(trimeshSolid) 
            sumVertexDefects, maxVertexDefects, maxAngleBetweenNormals = calculateNormalsDiff(trimeshSolid)  
            normalsPenalty = c2 * sumVertexDefects
            maxNormalsTerm = c3 * maxAngleBetweenNormals
            maxVertexDefectsTerm = c4 * maxVertexDefects    

            # print(f"Sum of vertex defects: {sumVertexDefects}")

            return [maxHFTerm + sumHFTerm + normalsPenalty + maxNormalsTerm + maxVertexDefectsTerm, normalsPenalty, maxNormalsTerm]    
            

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
                directoryName = self.makeDirectories(f"newspheretest{idx_to_vary}/1/test_", coefficients_list)
                #meshHFOpt(self, hfObjectiveFcn, constraint, updateHFProfile, calcHFAllMesh, calcMaxHF, calcEnergy, meshObj, coefficientsList, threshold, delta, id):
                self.opt.meshHFOpt(
                    objectiveFunction,  
                    facesToKeep,
                    self.fwd.makeHFProfile,
                    self.fwd.calculateAllHF,
                    self.fwd.filteredCalculateMaxHF, #self.fwd.calculateMaxHF,
                    self.fwd.calculateIntegratedEnergy,
                    my_trimeshSolid, 
                    coefficients_list,
                    threshold=0.00000001, 
                    delta=0.01, 
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
                        objectiveFunction,  
                        facesToKeep,
                        self.fwd.makeHFProfile,
                        self.fwd.calculateAllHF,
                        self.fwd.filteredCalculateMaxHF, #self.fwd.calculateMaxHF,
                        self.fwd.calculateIntegratedEnergy,
                        my_trimeshSolid, 
                        coefficients_list,
                        threshold=0.00000001, 
                        delta=0.01, 
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

        ## For bulk variable sweep testing
        # sweep_values_c0 = [0, 50, 100, 150]
        # sweep_values_c1 = [0, 50, 100, 150]
        # sweep_values_c2 = [0, 50, 100, 150]
        # sweep_values_c3 = [0, 50, 100, 150]
        # big_sweep_coefficients_and_record_output(sweep_values_c0, sweep_values_c1, sweep_values_c2, sweep_values_c3, id=0)

        # sweep_values_c0 = [100, 200, 300, 400]
        # sweep_values_c1 = [100, 200, 300, 400]
        # sweep_values_c2 = [100, 200, 300, 400]
        # sweep_values_c3 = [100, 200, 300, 400]
        # big_sweep_coefficients_and_record_output(sweep_values_c0, sweep_values_c1, sweep_values_c2, sweep_values_c3, id=1)

        # sweep_values_c0 = [200, 350, 500]
        # sweep_values_c1 = [350, 400, 450]
        # sweep_values_c2 = [500, 1000, 1500]
        # sweep_values_c3 = [500, 2000, 4000]
        # big_sweep_coefficients_and_record_output(sweep_values_c0, sweep_values_c1, sweep_values_c2, sweep_values_c3, id=3)


        ## For variable sweep testing
        # coefficients_list = [0, 0, 1500, 0, 0]
        # coefficients_list = [0, 0, 29698.489, 0, 0]
        # coefficients_list = [0, 0, 29698.48962 * 84.8528, 0, 0]
        # coefficients_list = [0, 0, 29698.48962, 0, 0] #original was * normalsPenalty_initial
        # c2val = [50, 100, 500, 1000] #(5000 * 504)/84.85281374238572
        # coefficients_list = [0, 0, c2val, 0, 0]

        # sweep_c3 = [200]
        #overall : c2 = c2 * sumVertexDefects, c3 = c3 * maxAngleBetweenNormals, c4 = c4 * maxVertexDefects   
        # coefficients_list = [0, 0, 0, 0, 0]
        # # sweep_c2 = [50, 500, 5000]
        # sweep_c2 = [46, 48, 52]
        # sweep_coefficients_and_record_output(coefficients_list, 2, sweep_c2)
        # coefficients_list = [0, 0, 50, 0, 0]
        # sweep_c3 = [10, 20, 30, 40]
        # sweep_coefficients_and_record_output(coefficients_list, 3, sweep_c3)

        # coefficients_list = [0, 0, 50, 0, 0]
        # sweep_c4 = [10, 50, 100, 500]
        # sweep_coefficients_and_record_output(coefficients_list, 4, sweep_c4)    

        coefficients_list = [0, 0, 50, 0, 10]
        sweep_c3 = [10, 50, 100, 500]
        sweep_coefficients_and_record_output(coefficients_list, 3, sweep_c3)   

        # coefficients_list = [0, 0, 0, 200, 0]
        # # sweep_c2 = [1000, 5000, 10000]
        # import random 
        # sweep_c2 = [random.uniform(800, 10000) for _ in range(3)]
        # sweep_coefficients_and_record_output(coefficients_list, 2, sweep_c2)

        # import random
        # sweep_values_c0 = [0]
        # sweep_values_c1 = [0]
        # # sweep_values_c2 = [random.uniform(700, 7000) for _ in range(2)]
        # sweep_values_c2 = [random.uniform(10000, 34000) for _ in range(2)]
        # sweep_values_c3 = [random.uniform(100, 800) for _ in range(2)]
        # big_sweep_coefficients_and_record_output(sweep_values_c0, sweep_values_c1, sweep_values_c2, sweep_values_c3, id=0)

        return 


if __name__ == '__main__':

    t0 = time.time()

    setup = RunSetup_MeshHF()

    """
    Run cube mesh operation test
    """
    setup.runOptimization()

    # p = Process(target = setup.runOptimization, args=())
    # p.start()
    # p.join()

    print(f"Time elapsed: {time.time() - t0}")





        
