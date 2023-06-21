import numpy as np

class ForwardModel_MeshHF: 

    def __init__(self, solidObj, q_mag, q_dir): 
        self.solidObj = solidObj
        self.q_mag = q_mag #magnitude of applied q [W/m^2]
        self.q_dir = q_dir #direction of applied q, [x,y,z] [m]
        return
    
    def calculateHFMeshElements(self, trimeshSolid):
        """
        Calculate heat flux on every mesh element and return them all in a list
        """
        normals, centers, areas = self.solidObj.normalsCentersAreas_Trimesh(trimeshSolid)
        # print(f"Normal vectors: {normals}")
        q_mesh_all = []
        for i in range(len(normals)):
            n = normals[i] 
            # print(f"Normal vector: {n}")
            q_i = np.dot(self.q_dir, n) * self.q_mag
            q_mesh_all.append(q_i)
        return q_mesh_all
    
    def calculateHFMeshSum(self, trimeshSolid):
        """
        Calculate sum of heat flux from all mesh elements
        """
        q_mesh_sum = 0
        normals, centers, areas = self.solidObj.normalsCentersAreas_Trimesh(trimeshSolid)

        for i in range(len(normals)):
            
            n = normals[i] 
            #dotprod >= 0.0 and dotprod <= 1.0
            dotprod = np.dot(self.q_dir, n)
            if dotprod > 0.0 and dotprod <= 1.0: #back face culling
                q_i = dotprod * self.q_mag
                q_mesh_sum += q_i

        #return np.sum(self.calculateHFMeshElements(trimeshSolid))
        return q_mesh_sum
    
    def calculateMaxHF(self, trimeshSolid):
        """
        Find single highest heat flux from all mesh element heat fluxes
        """
        q_mesh_all = self.calculateHFMeshElements(trimeshSolid)
        maxHF = np.max(q_mesh_all)
        return maxHF
    

from pymoo.optimize import minimize
from pymoo.model.problem import Problem
from pymoo.algorithms.so_genetic_algorithm import GA
from pymoo.factory import get_sampling, get_crossover, get_mutation
from pymoo.factory import get_termination
import trimesh

"""
for approach where i attempt using GA
""" 
class ForwardModelGA(Problem):

    def __init__(self, originalTriMesh, calculateHF):
        super().__init__(n_var=originalTriMesh.vertices.size, 
                         n_obj=1, 
                         xl=-np.inf, 
                         xu=np.inf)
        self.mesh = originalTriMesh
        self.calculateHF = calculateHF
        print("forward model set up")
        return
    
    def _evaluate(self, x, out, *args, **kwargs):
        """
        objective function - this is fitness criteria for GA
        """
        new_mesh = trimesh.Trimesh(vertices=x.reshape(self.mesh.vertices.shape), faces=self.mesh.faces)
        flux = self.calculateHF(new_mesh)
        print(f"flux calculated: {flux}")
        out["F"] = flux

