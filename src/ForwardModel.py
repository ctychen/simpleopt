import numpy as np

class ForwardModel_MeshHF: 

    def __init__(self, solidObj, q_mag, q_dir): 
        self.solidObj = solidObj
        self.q_mag = q_mag #magnitude of applied q [W/m^2]
        self.q_dir = q_dir #direction of applied q, [x,y,z] [m]
        return
    
    def calculateHFMesh(self, trimeshSolid):
        normals, centers, areas = self.solidObj.normalsCentersAreas_Trimesh(trimeshSolid)
        print(f"Normal vectors: {normals}")
        q_mesh_all = []
        for i in range(len(normals)):
            n = normals[i] 
            print(f"Normal vector: {n}")
            q_i = np.dot(self.q_dir, n) * self.q_mag
            q_mesh_all.append(q_i)
        return q_mesh_all
    
    def calculateMaxHF(self, trimeshSolid):
        q_mesh_all = self.calculateHFMesh(trimeshSolid)
        maxHF = np.max(q_mesh_all)
        return maxHF
    
    