import time
import sys
import numpy as np
import os

import pandas as pd

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

import matplotlib.pyplot as plt

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

        stlPath = " " #"box.stl"
        qDirIn = [0.0, -1.0, 0.0] #[m]
        # qDirIn = [0.707, -0.707, 0.0] #[m]
        # qDirIn = [-0.707, -0.707, 0.0] #[m]
        # qDirIn = [0.0, -0.707, 0.707] #[m]
        # qDirIn = [0.0, -0.707, -0.707] #[m]
        qMagIn = 10.0 #[W/m^2]

        self.box = Solid.MeshSolid(stlPath, stpPath)
        #self.box = Solid.MeshSolid(stlPath, stpPath) #normally, use this one!

        self.fwd = ForwardModel.ForwardModel_MeshHF(self.box, qMagIn, qDirIn, hfMode='uniform') 
        # self.fwd = ForwardModel.ForwardModel_MeshHF(self.box, qMagIn, qDirIn, hfMode='exponnorm') 
        self.opt = OptModel.OptModel_MeshHF()

        return
    

    def makeDirectories(self, runID, coefficientsList):
        c1 = coefficientsList[0]
        c2 = coefficientsList[1]
        c3 = coefficientsList[2]
        c4 = coefficientsList[3]
        c5 = coefficientsList[4]
        runName = runID + f'_c1_{c1:.2f}_c2_{c2:.2f}_c3_{c3:.2f}_c4_{c4:.2f}_c5_{c5:.2f}' #f'_c1_{c1:.2f}_c2_{c2:.2f}_c3_{c3:.2f}_c4_{c4:.2f}'
        runName = runName.replace(".", "-")

        directoryName = f"{runName}" #running this within docker container means can't save to external without bindmount aaa

        os.makedirs(directoryName)
        # os.makedirs(f"{directoryName}/images")

        print(f"Made directory: {directoryName}")

        #where we are saving generated vtk's
        # directoryName = f"{directoryName}/vtks"
        # os.makedirs(directoryName)

        return directoryName

    def makeSweepCSV(self, c1List, c2List, c3List, c4List, maxhfList, runID):
        data = {
            'c1': c1List,
            'c2': c2List,
            'c3': c3List,
            'c4': c4List,
            'maxhf': maxhfList
        }
        df = pd.DataFrame(data)
        df.to_csv(f"{runID}_results.csv", index=None, header=True)

        return 

    def runOptimization(self):
        
        def calculateNormalsDiff(trimeshSolid):
            """
            calculate differences between normal vectors of adjacent faces - ideally this should also be minimized?  
            """
            #"(n, 2) list of face indices. Each pair of faces in the list shares an edge, making them adjacent"
            adjacency = trimeshSolid.face_adjacency 

            normals = trimeshSolid.face_normals
            normalsDiff = normals[adjacency[:, 0]] - normals[adjacency[:, 1]]

            normalsDiffMagnitude = np.linalg.norm(normalsDiff, axis=1)

            reference_direction = np.array([0, 1, 0])  #"upwards" direction that we want it to move in 
            normalRefDotProducts = np.dot(normals, reference_direction)

            return normalsDiffMagnitude, normalRefDotProducts


        self.box.processSolid()
        trimeshSolid = self.box.trimeshSolid
        # trimeshSolid.export(f"{directoryName}/initial.stl")

        def findConstrainedFaces(mesh_center_xvals, mesh_center_yvals, mesh_center_zvals):
            return np.where((mesh_center_yvals == 10.0))[0]
            #return np.where((mesh_center_yvals == 10.0) | (mesh_center_xvals == 0.0) | (mesh_center_xvals == 10.0) | (mesh_center_zvals == 0.0) | (mesh_center_zvals == 10.0))[0]
        
        def calculateHeatFluxDiff(trimeshSolid):
            adjacency_info = trimeshSolid.face_adjacency

            allHF = self.fwd.calculateAllHF(trimeshSolid)
            heat_flux_diffs = []

            for face_pair in adjacency_info:
                face1, face2 = face_pair

                # Get the heat flux values for the adjacent faces
                heat_flux_face1 = allHF[face1]
                heat_flux_face2 = allHF[face2]

                heat_flux_diff = abs(heat_flux_face1 - heat_flux_face2)
                heat_flux_diffs.append(heat_flux_diff)
            
            return heat_flux_diffs

        
        def objectiveFunction(trimeshSolid, coefficientsList, unconstrainedFaces):

            c1 = coefficientsList[0] #max heat flux term
            c2 = coefficientsList[1] #sum heat flux, unweighted, term
            c3 = coefficientsList[2] #normals diff term
            c4 = coefficientsList[3] #energy term

            c5 = coefficientsList[4] #heat flux diff term

            q_mesh_all = self.fwd.calculateAllHF(trimeshSolid)

            # maxHFTerm = c1 * self.fwd.filteredCalculateMaxHF(q_mesh_all) #self.fwd.filteredCalculateMaxHF(trimeshSolid, unconstrainedFaces), 
            maxHFTerm = c1 * self.fwd.filteredCalculateMaxHF(q_mesh_all, unconstrainedFaces)
            sumHFTerm = c2 * self.fwd.calculateHFMeshSum(q_mesh_all) #self.fwd.calculateHFMeshSum(trimeshSolid)

            normalsDiff, normalRefDotProducts = calculateNormalsDiff(trimeshSolid)
            normalsPenalty = c3 * np.sum(normalsDiff)

            # energyTerm = c4 * self.fwd.calculateIntegratedEnergy(q_mesh_all) #self.fwd.calculateIntegratedEnergy(trimeshSolid)
            energyTerm = c4 * self.fwd.calculateIntegratedEnergy(q_mesh_all, trimeshSolid)

            hfDiffTerm = 0 #c5 * np.sum(calculateHeatFluxDiff(trimeshSolid))

            return maxHFTerm + sumHFTerm + normalsPenalty + energyTerm + hfDiffTerm        
            

        def sweep_coefficients_and_record_output(coefficients_list, idx_to_vary, sweep_values):
            """
            run through possible values of one coefficient while keeping the rest fixed
            do optimization and save results and attempted combinations to csv
            """

            c1_runvals = []
            c2_runvals = []
            c3_runvals = []
            c4_runvals = []
            maxhf_vals = []

            for val in sweep_values:
                my_trimeshSolid = trimesh.copy()
                coefficients_list[idx_to_vary] = val
                directoryName = self.makeDirectories(f"sweep_c{idx_to_vary}", coefficients_list)
                maxHF = self.opt.meshHFOpt(
                    objectiveFunction,  
                    findConstrainedFaces,
                    self.fwd.calculateAllHF,
                    self.fwd.filteredCalculateMaxHF, #self.fwd.calculateMaxHF,
                    self.fwd.calculateIntegratedEnergy,
                    my_trimeshSolid, 
                    coefficientsList,
                    threshold=0.000001, 
                    delta=0.01, 
                    id=directoryName
                )[0]
                c1_runvals.append(coefficients_list[0])
                c2_runvals.append(coefficients_list[1])
                c3_runvals.append(coefficients_list[2])
                c4_runvals.append(coefficients_list[3])
                maxhf_vals.append(maxHF)
            
            self.makeSweepCSV(c1_runvals, c2_runvals, c3_runvals, c4_runvals, maxhf_vals, f"sweep_c{idx_to_vary}")
            return 


        # #more runs 
        my_trimeshSolid = trimeshSolid.copy()
        # coefficientsList = [21.16, 0.53, 8.95, c4]
        # coefficientsList = [0, 0, 0, c4]
        # coefficientsList = [30.0, 1.0, 14.0, 4.55]
        # coefficientsList = [21.16, 0.53, 14.0, 4.55]
        coefficientsList = [21.16, 0.53, 14.0, 4.55, 0.0]
        # directoryName = self.makeDirectories(f"sweep_c5_{self.fwd.q_dir[0]}_{self.fwd.q_dir[1]}_{self.fwd.q_dir[2]}", coefficientsList)
        directoryName = self.makeDirectories(f"profiletest0", coefficientsList)
        maxHF = self.opt.meshHFOpt(
                objectiveFunction,
                findConstrainedFaces,  
                self.fwd.calculateAllHF,
                self.fwd.filteredCalculateMaxHF, #self.fwd.calculateMaxHF,
                self.fwd.calculateIntegratedEnergy,
                my_trimeshSolid, 
                coefficientsList,
                threshold=0.000001, 
                delta=0.01, 
                id=directoryName
        )[0]
        
        return maxHF


if __name__ == '__main__':

    t0 = time.time()

    setup = RunSetup_MeshHF()

    """
    Run cube mesh operation test
    """
    setup.runOptimization()

    print(f"Time elapsed: {time.time() - t0}")





        



