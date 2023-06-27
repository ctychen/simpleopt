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

    def runOptimization(self, runID="0011_03"):

        os.makedirs(f"test{runID}")

        self.box.processSolid()
        trimeshSolid = self.box.trimeshSolid
        trimeshSolid.export(f"test{runID}/initial.stl")

        def constraintFcn(trimeshSolid):
            #return true if maxHF is below threshold, false if it's above threshold. 
            #we want to keep running opt until we reach solns that satisfy this
            return (self.fwd.calculateMaxHF(trimeshSolid) < 9.0)
        
        def newObjective(trimeshSolid):
            #giving this an attempt since the distribution is really shifted
            return self.fwd.meanHF(trimeshSolid) + self.fwd.stdHF(trimeshSolid)
        
        def compoundObjective(trimeshSolid):
            #weighting for now, lets try it: 
            #max val of maxHF is on the order of 10
            #max val of sumHF is on the order of 1000 so scale down to same range
            c1 = 0.6
            c2 = 0.4
            return c1*self.fwd.calculateMaxHF(trimeshSolid) + c2*self.fwd.calculateHFMeshSum(trimeshSolid)
        
        def calculateNormalsDiff(trimeshSolid):
            """
            calculate differences between normal vectors of adjacent faces - ideally this should also be minimized?  
            """
            #"(n, 2) list of face indices. Each pair of faces in the list shares an edge, making them adjacent"
            adjacency = trimeshSolid.face_adjacency 

            normals = trimeshSolid.face_normals
            normalsDiff = normals[adjacency[:, 0]] - normals[adjacency[:, 1]]

            # num_adjacentFaces = np.bincount(adjacency.flatten())
            # internalFaces = np.where(num_adjacentFaces == 3)[0]
            # internalFacesOnly = np.isin(adjacency, internalFaces).all(axis=1)

            # normalsDiff = normalsDiff[internalFacesOnly]
            normalsDiffMagnitude = np.linalg.norm(normalsDiff, axis=1)

            return normalsDiffMagnitude
        

        def objectiveFunction(trimeshSolid):
            c1 = 0.6 #0.0 #0.6
            c2 = 0.4 #0.0 #0.4
            maxHFTerm = c1*self.fwd.calculateMaxHF(trimeshSolid)
            sumHFTerm = c2*self.fwd.calculateHFMeshSum(trimeshSolid)

            c3 = 0.5 #1.0 #0.5
            normalsDiff = calculateNormalsDiff(trimeshSolid)
            normalsPenalty = np.sum(normalsDiff) * c3
            
            return maxHFTerm + sumHFTerm + normalsPenalty
        

# def objectiveFunctionWithNormalsConstraint(tri_mesh, smooth_surface, max_normals_change):
#     # Calculate some measure of the heat flux (your original objective function)
#     heat_flux = ...

#     # Calculate the distance of each point to the smooth surface
#     distances = np.linalg.norm(tri_mesh.vertices - smooth_surface.evaluate(tri_mesh.vertices), axis=1)

#     # Calculate the penalty term
#     penalty = np.sum(distances**2)

#     # Calculate the normals change
#     normals_change = calculate_normals_change(tri_mesh)

#     # Add a penalty if the max normals change exceeds the allowed max
#     if np.max(normals_change) > max_normals_change:
#         penalty += np.sum(normals_change)

#     # Return the objective function value
#     return heat_flux + lambda * penalty


        #optimizer setup
        # return self.opt.meshHFOpt(self.fwd.calculateHFMeshSum, trimeshSolid, self.opt.moveMeshVertices, threshold=100, delta=0.05, id=runID)
        
        #args:
        #hfObjectiveFcn, calcHFAllMesh, calcMaxHF, calcHFSum, meshObj, changeMeshFcn, threshold, delta, id

        return self.opt.meshHFOpt(
            objectiveFunction, #compoundObjective, #self.fwd.calculateHFMeshSum, #compoundObjective, 
            self.fwd.calculateAllHF,
            self.fwd.calculateMaxHF, #self.fwd.calculateMaxHF, #using this for now to plot components of compound obj
            self.fwd.calculateHFMeshSum,
            trimeshSolid, 
            self.opt.moveMeshVertices, 
            threshold=0.000001, 
            delta=0.01, 
            id=runID
        )


if __name__ == '__main__':

    t0 = time.time()

    setup = RunSetup_MeshHF()

    """
    Run cube mesh operation test
    """
    setup.runOptimization()

    print(f"Time elapsed: {time.time() - t0}")





        



