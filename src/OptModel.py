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
    

    def gradientDescentHF(self, tri_mesh, objectiveFunction, delta, fileID):
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
                # print(f"After moving: {tri_mesh.vertices[i, j]}")
                obj_afterMoving = objectiveFunction(tri_mesh)
                # print(f"Objective function after moving: {obj_afterMoving}")
                tri_mesh.export(f"test{fileID}/wip/movingvertices_{i}_{j}.stl")
                # print(f"Saved file for solid after moving vertices {i}, {j}")
                #move vertex back to original, calc objective function then for comparison
                tri_mesh.vertices[i, j] -= delta
                obj_beforeMoving = objectiveFunction(tri_mesh)
                # OR: maybe the change is to change the objective function so it's only calculating the HF on this one mesh element?
                #and then modifying one element at a time?
                # print(f"Before moving: {tri_mesh.vertices[i, j]}")
                # print(f"Objective function before moving: {obj_beforeMoving}")
                gradient[i, j] = (obj_afterMoving - obj_beforeMoving) / (2 * delta)
                # print(f"Gradient value: {gradient[i, j]}")
                # print(f"Gradient should be: {(obj_afterMoving - obj_beforeMoving) / (2 * delta)}")
        
        return gradient
    
    #self.optHFandModifyMesh(self, trimeshSolid, hfFunction, hfObjectiveFcn, stepSize=step, id=id, count=count, maxMargin=0.05)
    def optHFandModifyMesh(self, tri_mesh, objectiveFunction, stepSize, id, count):

        #q_mesh_all = hfFunction(tri_mesh) #fwdModel.calculateHFMeshElements(tri_mesh)

        #max_faces = np.where(q_mesh_all >= (np.max(q_mesh_all) - maxMargin))[0]

        #for each face that contributes to the maximum, calculate the gradient
        #of the objective function with respect to the positions of its vertices
        for face in tri_mesh.faces:
            # Get the vertices of this face
            vertices = tri_mesh.faces[face]

            #calculate the gradient for each vertex
            for vertex in vertices:
                
                grad = self.compute_gradient_single_vertex(tri_mesh, objectiveFunction, vertex)
                
                #update the position of this vertex
                #mesh.vertices[vertex] -= step_size * grad
                tri_mesh.vertices[vertex] -= stepSize * grad

                if count % 2 == 0: 
                    tri_mesh.export(f"test{id}/wip/{count}_hf_{objectiveFunction(tri_mesh)}.stl")

        return 
    

    def compute_gradient_single_vertex(self, mesh, objective_function, vertex_index, epsilon=0.5):
        # Store the original position
        original_pos = mesh.vertices[vertex_index].copy()

        # Initialize gradient
        gradient = np.zeros(3)

        origMesh = mesh

        # For each dimension
        for i in range(3):
            #move vertex slightly in this dimension
            mesh.vertices[vertex_index][i] += epsilon

            #calculate change in the objective function
            delta = objective_function(mesh) - objective_function(origMesh)

            #estimate gradient in this dimension
            gradient[i] = delta / epsilon

            #reset vertex to its original position to repeat process
            mesh.vertices[vertex_index] = original_pos

        return gradient

    
    def moveMeshVertices(self, trimeshSolid, gradient, delta):
        """
        function for how we want to adjust mesh vertices, depending on what the gradient is 
        """
        # return trimeshSolid.vertices - delta
        return trimeshSolid.vertices - (delta * gradient)

    def meshHFOpt(self, hfObjectiveFcn, meshObj, changeMeshFcn, threshold, delta, id):
    # def meshHFOpt(self, hfFunction, hfObjectiveFcn, meshObj, threshold, step, id):
        """
        runs optimization process until objective fcn value reaches stopping condition @ minimum
        modifies the mesh based on gradient by applying changeMeshFcn accordingly

        can change changeMeshFcn, hfObjectiveFcn to use different functions
        if we want a different manipulation, or add more stuff to the functions
        """
        #assuming input is already a trimesh, ie. processing the solid was done already
        trimeshSolid = meshObj
        count = 0

        print("Starting the mesh HF opt")

        print(f"Starting max HF: {hfObjectiveFcn(trimeshSolid)}")

        #objective fcn should take in a trimesh mesh object
        while hfObjectiveFcn(trimeshSolid) > threshold:

            print(f"Current max HF: {hfObjectiveFcn(trimeshSolid)}")

            #calc the gradient
            gradient = self.gradientDescentHF(trimeshSolid, hfObjectiveFcn, delta, fileID=id)
            print(f"Gradient calculated: {gradient}")

            #move the vertices a bit based on the gradient
            trimeshSolid.vertices = changeMeshFcn(trimeshSolid, gradient, delta)

            trimeshSolid.export(f"test{id}/wip_{count}.stl")

            #optHFandModifyMesh(self, tri_mesh, hfFunction, objectiveFunction, stepSize, id, count, maxMargin=0.05):
            #def optHFandModifyMesh(self, tri_mesh, objectiveFunction, stepSize, id, count):
            # self.optHFandModifyMesh(trimeshSolid, hfObjectiveFcn, stepSize=step, id=id, count=count)
            #self.optHFandModifyMesh(trimeshSolid, hfFunction, hfObjectiveFcn, stepSize=step, id=id, count=count, maxMargin=0.5)
            print(f"Current max HF: {hfObjectiveFcn(trimeshSolid)}")
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

