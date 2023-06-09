import numpy as np

class ForwardModel_Box: 
    def __init__(self, g_x, solid, q_mag, q_dir): #some arbitrary starting values todo
        self.solid = solid
        self.q_mag = q_mag #magnitude of applied q [W/m^2]
        self.q_dir = q_dir #direction of applied q, [x,y,z] [m]
        self.g_x = g_x #objective function, can swap out
        return

    #right now it's simple but ultimately this fcn could have whatever in it - formerly processCADModel()
    def process(self):
        self.solid.processModel()
        return
    
    def calcQMesh(self):
        #for every mesh element, calculate q and store in list/arr with same length as mesh elem list
        q_mesh_all = []
        for i in range(len(self.norms[0])): 
            norm = self.norms[0][i]
            # print(norm)
            q_i = np.dot(self.q_dir, norm) * self.q_mag
            q_mesh_all.append(q_i)
        # print(q_mesh_all)
        return q_mesh_all
    
    def calcQMesh_Vector(self): #version to use for list of vertices
        q_mesh_all = []
        for i in range(len(self.solid.norms)): 
            norm = self.solid.norms[i]
            # print(norm)
            q_i = np.dot(self.q_dir, norm) * self.q_mag
            q_mesh_all.append(q_i)
        return q_mesh_all
    
    def calcQMesh_Vector(self, vertices):
        q_mesh_all = []
        norms = self.solid.normsCentersAreas_VectorAny(vertices)[0]
        #want to calculate norms of the generated mesh... 

        for norm in norms:
            dotprod = np.dot(self.q_dir, norm)
            if dotprod >= 0.0 and dotprod <= 1.0:
                q_i = dotprod * self.q_mag
                q_mesh_all.append(q_i) 

        # print(f"Q mesh all: {q_mesh_all}")
        return q_mesh_all
    
    def calcQMesh_Vec(self, vector):
        return np.dot(self.fwd.q_dir, vector) * self.fwd.q_mag
    
    def calcObjective(self, q_mesh_all): 
        #takes in objective function
        return self.g_x(q_mesh_all)
    
from pymoo.model.problem import Problem

class ForwardModel_MeshOptimization(Problem):
    
    def __init__(self, vertices, solid, n_vars):
        super().__init__(n_var=n_vars, n_obj=1, n_constr=0, xl=-1, xu=1)  # Set variable bounds and other details
        self.vertices = vertices #maybe don't set this here, and get from property of Solid???
        self.solid = solid

    """
    Evaluate HF on entire modified mesh
    """
    def _evaluate(self, x, out, *args, **kwargs):
        #this should return a set of solutions though - min(q_mesh_all) would be just 1?
        #calculate HF on entire modified mesh - eval quality
        modifiedVertices = self.solid.modifyMesh(self.vertices, x)
        q_mesh_all = self.calculateQMeshAll(modifiedVertices)
        out["F"] = q_mesh_all

    """
    q_mag: [W/m^2]
    q_dir: [m], (x,y,z)
    """
    def calculateQMeshAll(self, vertices, q_mag, q_dir):
        q_mesh_all = []

        norms = self.solid.normsCentersAreas_VectorAny(vertices)[0]

        for norm in norms:
            dotprod = np.dot(q_dir, norm)
            if dotprod >= 0.0 and dotprod <= 1.0:
                q_i = dotprod * q_mag
                q_mesh_all.append(q_i) 
        
        return q_mesh_all