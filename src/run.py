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
        
        
        c1 = np.random.rand() * 10 #5.0 #for maxHF
        c2 = np.random.rand() / 1.5 #0.1 #for sum HF
        c3 = np.random.rand() #/ 2.0 #0.5 #for normals of surfaces
        c4 = np.random.rand() / 1.5 #0.2 #for energy
        c5 = np.random.rand() #for distances from original mesh

        runName = runID + f'_c1_{c1:.2f}_c2_{c2:.2f}_c3_{c3:.2f}_c4_{c4:.2f}_c5_{c5:.2f}'  #runID + f"_c1_{c1.2f}_c2_{c2:03}_c3_{c3:03}_c4_{c4:03}"
        runName = runName.replace(".", "-")

        directoryName = f"{runName}" #running this within docker container means can't save to external without bindmount aaa

        os.makedirs(directoryName)
        os.makedirs(f"{directoryName}/images")

        print(f"Made directory: {directoryName}")

        #where we are saving generated vtk's
        directoryName = f"{directoryName}/vtks"
        os.makedirs(directoryName)

        self.box.processSolid()
        trimeshSolid = self.box.trimeshSolid
        trimeshSolid.export(f"{directoryName}/initial.stl")

        originalTrimesh = trimeshSolid

        def objectiveFunction(trimeshSolid):
            # c1 = 5.0 #0.6
            # c2 = 0.1 #0.1 #1.0 #0.0 #0.4
            maxHFTerm = c1*self.fwd.calculateMaxHF(trimeshSolid)
            sumHFTerm = c2*self.fwd.calculateHFMeshSum(trimeshSolid)

            # c3 = 0.5 #0.0 #1.0 #0.5
            normalsDiff = calculateNormalsDiff(trimeshSolid)
            normalsPenalty = np.sum(normalsDiff) * c3

            # c4 = 0.2
            energyTerm = c4*self.fwd.calculateIntegratedEnergy(trimeshSolid)

            #calc distance of each vertex to the original mesh
            #calc penalty for being too close to the original surface ()
            #trimesh signed distances: 
            #points OUTSIDE the mesh will have NEGATIVE distance
            #points within tol.merge of the surface will have POSITIVE distance
            #points INSIDE the mesh will have POSITIVE distance
            #and overall we want to force the mesh to move to OUTSIDE the original mesh if we don't want it punching through bottom
            distances = trimesh.proximity.signed_distance(originalTrimesh, trimeshSolid.vertices)
            distancePenalty = -np.sum(np.abs(distances))
            distancesTerm = c5*distancePenalty
            
            return maxHFTerm + sumHFTerm + normalsPenalty + energyTerm + distancesTerm

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
            id=directoryName#runName #runID
        )


if __name__ == '__main__':

    t0 = time.time()

    setup = RunSetup_MeshHF()

    """
    Run cube mesh operation test
    """
    setup.runOptimization()

    print(f"Time elapsed: {time.time() - t0}")





        



