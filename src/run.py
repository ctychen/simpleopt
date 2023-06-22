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

        stpPath = "unit_test_cone.step" #"unit_test_cube.step"

        stlPath = " " #"box.stl"
        qDirIn = [0.0, -1.0, 0.0] #[m]
        qMagIn = 10.0 #[W/m^2]

        self.box = Solid.MeshSolid(stlPath, stpPath)
        #self.box = Solid.MeshSolid(stlPath, stpPath) #normally, use this one!

        self.fwd = ForwardModel.ForwardModel_MeshHF(self.box, qMagIn, qDirIn) 
        self.opt = OptModel.OptModel_MeshHF()

        return

    def runOptimization(self, runID="006_3"):

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

        #optimizer setup
        # return self.opt.meshHFOpt(self.fwd.calculateHFMeshSum, trimeshSolid, self.opt.moveMeshVertices, threshold=100, delta=0.05, id=runID)
        
        return self.opt.meshHFOpt(
            #self.fwd.meanHF, 
            #newObjective,
            self.fwd.calculateHFMeshSum,
            constraintFcn, 
            self.fwd.calculateHFMeshElements, 
            #self.fwd.distForObj, 
            self.fwd.calculateMaxHF,
            trimeshSolid, 
            self.opt.moveMeshVertices, 
            threshold=0.8, 
            delta=0.05, 
            id=runID
            )

        # args: hfObjectiveFcn, meshObj, changeMeshFcn, threshold, delta
        # args: meshHFOpt(self, hfFunction, hfObjectiveFcn, meshObj, threshold, stepSize, id)
        #calculateHFMeshSum
        # def meshHFOpt(self, hfObjectiveFcn, meshObj, changeMeshFcn, threshold, delta, id):
        #return self.opt.meshHFOpt(self.fwd.calculateHFMeshElements, self.fwd.calculateMaxHF, trimeshSolid, threshold=0.1, step=0.1, id=runID)
        # return self.opt.meshHFOpt(self.fwd.calculateHFMeshElements, self.fwd.calculateHFMeshSum, trimeshSolid, threshold=0.1, step=0.1, id=runID)


if __name__ == '__main__':

    t0 = time.time()

    setup = RunSetup_MeshHF()

    """
    Run cube mesh operation test
    """
    setup.runOptimization()

    print(f"Time elapsed: {time.time() - t0}")





        



