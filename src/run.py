import time
import sys
import numpy as np
import os

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
        qMagIn = 10.0 #[W/m^2]

        self.box = Solid.MeshSolid(stlPath, stpPath)
        #self.box = Solid.MeshSolid(stlPath, stpPath) #normally, use this one!

        self.fwd = ForwardModel.ForwardModel_MeshHF(self.box, qMagIn, qDirIn) 
        self.opt = OptModel.OptModel_MeshHF()

        return

    def runOptimization(self, runID="0012"):

        os.makedirs(f"test{runID}")

        self.box.processSolid()
        trimeshSolid = self.box.trimeshSolid
        trimeshSolid.export(f"test{runID}/initial.stl")

        
        def calculateNormalsDiff(trimeshSolid):
            """
            calculate differences between normal vectors of adjacent faces - ideally this should also be minimized?  
            """
            #"(n, 2) list of face indices. Each pair of faces in the list shares an edge, making them adjacent"
            adjacency = trimeshSolid.face_adjacency 

            normals = trimeshSolid.face_normals
            normalsDiff = normals[adjacency[:, 0]] - normals[adjacency[:, 1]]


            ##below lines are if we want to detect faces that belonged to cube edges
            ##and then we can exclude those when we use the normals for the objective fcn
            # num_adjacentFaces = np.bincount(adjacency.flatten())
            # internalFaces = np.where(num_adjacentFaces == 3)[0]
            # internalFacesOnly = np.isin(adjacency, internalFaces).all(axis=1)
            # normalsDiff = normalsDiff[internalFacesOnly]

            normalsDiffMagnitude = np.linalg.norm(normalsDiff, axis=1)

            return normalsDiffMagnitude
        
        #np.random.rand() returns [0, 1)
        c1 = 5.0 #np.random.rand() * 10
        c2 = 0.1 #np.random.rand() / 2.0
        c3 = 0.5 #np.random.rand() / 2.0
        c4 = 0.2 #np.random.rand() * 10

        runName = runID + f"_c1_{c1}_c2_{c2}_c3_{c3}_c4_{c4}"

        def objectiveFunction(trimeshSolid):
            # c1 = 5.0 #0.6
            # c2 = 0.1 #0.1 #1.0 #0.0 #0.4
            maxHFTerm = c1*self.fwd.calculateMaxHF(trimeshSolid)
            sumHFTerm = c2*self.fwd.calculateHFMeshSum(trimeshSolid)

            c3 = 0.5 #0.0 #1.0 #0.5
            normalsDiff = calculateNormalsDiff(trimeshSolid)
            normalsPenalty = np.sum(normalsDiff) * c3

            c4 = 0.2
            energyTerm = c4*self.fwd.calculateIntegratedEnergy(trimeshSolid)
            
            return maxHFTerm + sumHFTerm + normalsPenalty + energyTerm

        #optimizer setup
        # return self.opt.meshHFOpt(self.fwd.calculateHFMeshSum, trimeshSolid, self.opt.moveMeshVertices, threshold=100, delta=0.05, id=runID)
        
        #args:
        #hfObjectiveFcn, calcHFAllMesh, calcMaxHF, calcHFSum, meshObj, changeMeshFcn, threshold, delta, id

        return self.opt.meshHFOpt(
            objectiveFunction, #compoundObjective, #self.fwd.calculateHFMeshSum, #compoundObjective, 
            self.fwd.calculateAllHF,
            self.fwd.calculateMaxHF, #self.fwd.calculateMaxHF, #using this for now to plot components of compound obj
            #self.fwd.calculateHFMeshSum,
            self.fwd.calculateIntegratedEnergy,
            trimeshSolid, 
            self.opt.moveMeshVertices, 
            threshold=0.000001, 
            delta=0.01, 
            id=runName #runID
        )


if __name__ == '__main__':

    t0 = time.time()

    setup = RunSetup_MeshHF()

    """
    Run cube mesh operation test
    """
    setup.runOptimization()

    print(f"Time elapsed: {time.time() - t0}")





        



