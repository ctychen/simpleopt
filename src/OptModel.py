import FreeCAD
import Part
import Mesh
import MeshPart
from FreeCAD import Base
import numpy as np

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots


class OptModel_MeshHF: 
    """
    Mesh element manipulation to optimize heat flux - model for optimization
    """

    def __init__(self):
        """
        TODO: figure out what properties we actually need for optmodel, and what can just be fcn inputs
        """
        return
    

    def gradientDescentHF(self, tri_mesh, objectiveFunction, delta):
        """
        gradient descent implementation for heat flux minimization
        takes in trimesh object, applies objective function to it
        moves all mesh vertices by a small amount and finds the gradient when doing so 
        """
        gradient = np.zeros_like(tri_mesh.vertices)
        #for each vertex
        for i in range(len(tri_mesh.vertices)):
            #for each dimension
            for j in range(3):
                #move vertex a bit, see if the objective function decreases
                tri_mesh.vertices[i, j] += delta
                obj_afterMoving = objectiveFunction(tri_mesh)
                #move vertex back to original, calc objective function then for comparison
                tri_mesh.vertices[i, j] -= delta
                obj_beforeMoving = objectiveFunction(tri_mesh)
                gradient[i, j] = (obj_afterMoving - obj_beforeMoving) / (2 * delta)
        return gradient
    
    def moveMeshVertices(self, trimeshSolid, gradient, delta):
        """
        function for how we want to adjust mesh vertices, depending on what the gradient is 
        """
        return trimeshSolid.vertices - (delta * gradient)

    def meshHFOpt(self, hfObjectiveFcn, meshObj, changeMeshFcn, threshold, delta, id):
        """
        runs optimization process until objective fcn value reaches stopping condition @ minimum
        modifies the mesh based on gradient by applying changeMeshFcn accordingly

        can change changeMeshFcn, hfObjectiveFcn to use different functions
        if we want a different manipulation, or add more stuff to the functions
        """
        #assuming input is already a trimesh, ie. processing the solid was done already
        trimeshSolid = meshObj
        count = 0

        #objective fcn should take in a trimesh mesh object
        while hfObjectiveFcn(trimeshSolid) > threshold:

            #calc the gradient
            #gradientDescentHF(self, tri_mesh, objectiveFunction, delta)
            gradient = self.gradientDescentHF(trimeshSolid, hfObjectiveFcn, delta)
            print(f"Gradient calculated: {gradient}")
            # gradient = gradientDescentFcn()

            #move the vertices a bit based on the gradient (not sure if you can do this without looping)
            #trimeshSolid.Vertices -= 0.01 * gradient

            trimeshSolid.vertices = changeMeshFcn(trimeshSolid, gradient, delta)

            trimeshSolid.export(f"test{id}/wip_{count}.stl")
            count += 1
        
        #when process is done, the mesh should have been modified - so return it 
        return trimeshSolid



class OptModel_Template:
    """
    Template class for optimizer, change these functions to match problem we're dealing with 
    """
    def __init__(self, threshold = 0.1, gprev = 10000, gcurr = 5, delstep = 0.1, del_e = 1):
        #all placeholder values for now I have no idea what they should be yet - maybe these shouldn't be here and instead in opt? 
        self.threshold_err = threshold
        self.g_prev = gprev
        self.g_curr = gcurr
        self.delstep = delstep
        self.del_e = del_e
        return
    
    def calculateDelE(self): 
        self.del_e = self.g_prev - self.g_curr
        return self.del_e

    
    def doTransform(self, cadModel, x=0, y=0, z=1): 
        return
    

    def updategValues(self, g_new):
        self.g_prev = self.g_curr
        self.g_curr = g_new
        return    

