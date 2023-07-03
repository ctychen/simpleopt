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
        qMagIn = 50.0 #10.0 #[W/m^2]

        self.box = Solid.MeshSolid(stlPath, stpPath)
        #self.box = Solid.MeshSolid(stlPath, stpPath) #normally, use this one!

        self.fwd = ForwardModel.ForwardModel_MeshHF(self.box, qMagIn, qDirIn) 
        self.opt = OptModel.OptModel_MeshHF()

        #to make nonuniform, eich-like HF profile on top face
        self.fwd.makeHFProfile(self.box.trimeshSolid, directionVector=qDirIn)

        return
    

    def makeDirectories(self, runID, coefficientsList):
        c1 = coefficientsList[0]
        c2 = coefficientsList[1]
        c3 = coefficientsList[2]
        c4 = coefficientsList[3]
        runName = runID + f'_c1_{c1:.2f}_c2_{c2:.2f}_c3_{c3:.2f}_c4_{c4:.2f}'
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

        # c1 = 21.16
        # c2 = 0.0 
        # c3 = 8.95
        # c4 = 4.55
        # c5 = 0.0 

        # runName = runID + f'_c1_{c1:.2f}_c2_{c2:.2f}_c3_{c3:.2f}_c4_{c4:.2f}_c5_{c5:.2f}'
        # runName = runName.replace(".", "-")

        # directoryName = f"{runName}" #running this within docker container means can't save to external without bindmount aaa

        # os.makedirs(directoryName)
        # os.makedirs(f"{directoryName}/images")

        # print(f"Made directory: {directoryName}")

        # #where we are saving generated vtk's
        # directoryName = f"{directoryName}/vtks"
        # os.makedirs(directoryName)

        self.box.processSolid()
        trimeshSolid = self.box.trimeshSolid
        # trimeshSolid.export(f"{directoryName}/initial.stl")
        
        # def objectiveFunction(trimeshSolid, unconstrainedFaces):
        def objectiveFunction(trimeshSolid, coefficientsList):

            c1 = coefficientsList[0]
            c2 = coefficientsList[1]
            c3 = coefficientsList[2]
            c4 = coefficientsList[3]

            # print(f"Coefficients used: {coefficientsList}")
            
            maxHFTerm = c1 * self.fwd.calculateMaxHF(trimeshSolid)
            sumHFTerm = c2 * self.fwd.calculateHFMeshSum(trimeshSolid)

            normalsDiff, normalRefDotProducts = calculateNormalsDiff(trimeshSolid)
            normalsPenalty = c3 * np.sum(normalsDiff)

            energyTerm = c4 * self.fwd.calculateIntegratedEnergy(trimeshSolid)

            #attempting to use this to see if we can encourage the mesh to go more convex?
            distancesTerm = 0.0
            # distancesTerm = c5 * np.sum(np.where(normalRefDotProducts < 0, -normalRefDotProducts, 0))

            # print(f"Value of objective: {maxHFTerm + sumHFTerm + normalsPenalty + energyTerm + distancesTerm    }")

            return maxHFTerm + sumHFTerm + normalsPenalty + energyTerm + distancesTerm            
            

        # coefficientsList = [21.16, 0.53, 8.95, 4.55]
        coefficientsList = [21.16, 1.00, 8.95, 4.55]
        my_trimeshSolid = trimeshSolid.copy()

        directoryName = self.makeDirectories("015_nonuniform", coefficientsList)
        maxHF, new_trimeshSolid = self.opt.meshHFOpt(
            objectiveFunction,  
            self.fwd.calculateAllHF,
            self.fwd.calculateMaxHF,
            self.fwd.calculateIntegratedEnergy,
            my_trimeshSolid, 
            coefficientsList,
            threshold=0.000001, 
            delta=0.01, 
            id=directoryName
        )[0]

        return maxHF, new_trimeshSolid



        # c1_runvals = []
        # c2_runvals = []
        # c3_runvals = []
        # c4_runvals = []
        # maxhf_vals = []

        # C1: 21.16
        # C2: 0.53
        # C3: 8.95
        # C4: 4.55

        # c1_sweep = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0]
        # for c1 in c1_sweep:
        #     my_trimeshSolid = trimeshSolid.copy()
        #     coefficientsList = [c1, 0.53, 8.95, 4.55]
        #     directoryName = self.makeDirectories("sweep_c1", coefficientsList)
        #     maxHF = self.opt.meshHFOpt(
        #         objectiveFunction,  
        #         self.fwd.calculateAllHF,
        #         self.fwd.calculateMaxHF,
        #         self.fwd.calculateIntegratedEnergy,
        #         my_trimeshSolid, 
        #         coefficientsList,
        #         threshold=0.000001, 
        #         delta=0.01, 
        #         id=directoryName
        #     )[0]
        #     c1_runvals.append(c1)
        #     c2_runvals.append(coefficientsList[1])
        #     c3_runvals.append(coefficientsList[2])
        #     c4_runvals.append(coefficientsList[3])
        #     maxhf_vals.append(maxHF)

        # self.makeSweepCSV(c1_runvals, c2_runvals, c3_runvals, c4_runvals, maxhf_vals, "sweep_c1")

        # # c2_sweep = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]
        # c2_sweep = [0.6, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0]
        # for c2 in c2_sweep:
        #     my_trimeshSolid = trimeshSolid.copy()
        #     coefficientsList = [21.16, c2, 8.95, 4.55]
        #     directoryName = self.makeDirectories("sweep_c2", coefficientsList)
        #     maxHF = self.opt.meshHFOpt(
        #         objectiveFunction,  
        #         self.fwd.calculateAllHF,
        #         self.fwd.calculateMaxHF,
        #         self.fwd.calculateIntegratedEnergy,
        #         my_trimeshSolid, 
        #         coefficientsList,
        #         threshold=0.000001, 
        #         delta=0.01, 
        #         id=directoryName
        #     )[0]
        #     c1_runvals.append(coefficientsList[0])
        #     c2_runvals.append(c2)
        #     c3_runvals.append(coefficientsList[2])
        #     c4_runvals.append(coefficientsList[3])
        #     maxhf_vals.append(maxHF)

        # self.makeSweepCSV(c1_runvals, c2_runvals, c3_runvals, c4_runvals, maxhf_vals, "sweep_c2")

        # c3_sweep = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0]
        # for c3 in c3_sweep:
        #     my_trimeshSolid = trimeshSolid.copy()
        #     coefficientsList = [21.16, 0.53, c3, 4.55]
        #     directoryName = self.makeDirectories("sweep_c3", coefficientsList)
        #     maxHF = self.opt.meshHFOpt(
        #         objectiveFunction,  
        #         self.fwd.calculateAllHF,
        #         self.fwd.calculateMaxHF,
        #         self.fwd.calculateIntegratedEnergy,
        #         my_trimeshSolid, 
        #         coefficientsList,
        #         threshold=0.000001, 
        #         delta=0.01, 
        #         id=directoryName
        #     )[0]
        #     c1_runvals.append(coefficientsList[0])
        #     c2_runvals.append(coefficientsList[1])
        #     c3_runvals.append(c3)
        #     c4_runvals.append(coefficientsList[3])
        #     maxhf_vals.append(maxHF)

        # self.makeSweepCSV(c1_runvals, c2_runvals, c3_runvals, c4_runvals, maxhf_vals, "sweep_c3")

        # c4_sweep = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
        # for c4 in c4_sweep:
        #     my_trimeshSolid = trimeshSolid.copy()
        #     coefficientsList = [21.16, 0.53, 8.95, c4]
        #     directoryName = self.makeDirectories("sweep_c4", coefficientsList)
        #     maxHF = self.opt.meshHFOpt(
        #         objectiveFunction,  
        #         self.fwd.calculateAllHF,
        #         self.fwd.calculateMaxHF,
        #         self.fwd.calculateIntegratedEnergy,
        #         my_trimeshSolid, 
        #         coefficientsList,
        #         threshold=0.000001, 
        #         delta=0.01, 
        #         id=directoryName
        #     )[0]
        #     c1_runvals.append(coefficientsList[0])
        #     c2_runvals.append(coefficientsList[1])
        #     c3_runvals.append(coefficientsList[2])
        #     c4_runvals.append(c4)
        #     maxhf_vals.append(maxHF)

        # self.makeSweepCSV(c1_runvals, c2_runvals, c3_runvals, c4_runvals, maxhf_vals, "sweep_c4")
        
        # return
    
        # my_trimeshSolid = trimeshSolid.copy()

        # return self.opt.meshHFOpt(
        #     objectiveFunction,  
        #     self.fwd.calculateAllHF,
        #     self.fwd.calculateMaxHF,
        #     self.fwd.calculateIntegratedEnergy,
        #     my_trimeshSolid, 
        #     coefficientsList,
        #     threshold=0.000001, 
        #     delta=0.01, 
        #     id=directoryName
        # )

        # return self.opt.meshHFOpt(
        #     objectiveFunction,  
        #     self.fwd.calculateAllHF,
        #     self.fwd.calculateMaxHF,
        #     self.fwd.calculateIntegratedEnergy,
        #     trimeshSolid, 
        #     threshold=0.000001, 
        #     delta=0.01, 
        #     id=directoryName
        # )


if __name__ == '__main__':

    t0 = time.time()

    setup = RunSetup_MeshHF()

    """
    Run cube mesh operation test
    """
    setup.runOptimization()

    print(f"Time elapsed: {time.time() - t0}")





        



