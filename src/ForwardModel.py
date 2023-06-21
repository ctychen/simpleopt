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
    
    