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

    def runOptimization(self, runID="009_5"):

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

        #optimizer setup
        # return self.opt.meshHFOpt(self.fwd.calculateHFMeshSum, trimeshSolid, self.opt.moveMeshVertices, threshold=100, delta=0.05, id=runID)
        
        return self.opt.meshHFOpt(
            #self.fwd.meanHF, 
            #newObjective,
            #self.fwd.calculateHFMeshSum,
            compoundObjective, #self.fwd.calculateHFMeshSum, #compoundObjective, 
            self.fwd.calculateAllHF,
            #self.fwd.calculateHFMeshElements, 
            #self.fwd.calculateHFMeshSum, #using this for now to plot components of compound obj
            #self.fwd.distForObj, 
            self.fwd.calculateMostNegHF, #self.fwd.calculateMaxHF, #using this for now to plot components of compound obj
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





        



