import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

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
        q_mesh_all = []
        for i in range(len(normals)):
            n = normals[i] 
            dotprod = np.dot(self.q_dir, n)
            # if dotprod > 0.0 and dotprod <= 1.0:
            
            if dotprod < 0.0 and abs(dotprod) <= 1.0: #not a backface/shadowed since shadowed would be if b dot n > 0
                q_i = abs(dotprod) * self.q_mag
                q_mesh_all.append(q_i)
        return q_mesh_all
    
    
    def calculateHFDistribution(self, trimeshSolid):
        #q_mesh_all = np.array(self.calculateHFMeshElements(trimeshSolid))
        normals, centers, areas = self.solidObj.normalsCentersAreas_Trimesh(trimeshSolid)
        q_mesh_all = []
        for i in range(len(normals)):
            n = normals[i] 
            dotprod = np.dot(self.q_dir, n)
            #adjusting the range here to try and get it to move not the edges by a bit - distribution is very much around 0 right now
            #if dotprod > 0.0 and dotprod <= 1.0: 

            if dotprod < 0.0 and abs(dotprod) <= 1.0:
                q_i = abs(dotprod) * self.q_mag
                q_mesh_all.append(q_i)

        q_mesh_all = np.array(q_mesh_all)

        hfMean = np.mean(q_mesh_all)
        hfVariance = np.var(q_mesh_all)
        hfStd = np.std(q_mesh_all)

        return hfMean, hfVariance, hfStd, q_mesh_all
    
    def meanHF(self, trimeshSolid):
        return self.calculateHFDistribution(trimeshSolid)[0]
    
    def varianceHF(self, trimeshSolid):
        return self.calculateHFDistribution(trimeshSolid)[1]
    
    def stdHF(self, trimeshSolid):
        return self.calculateHFDistribution(trimeshSolid)[2]
    
    def distForObj(self, trimeshSolid):
        return self.calculateHFDistribution(trimeshSolid)[3]
    
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
            #if dotprod > 0.0 and dotprod <= 1.0: 

            #this condition is bc highest fluxes would be expected to be on faces facing into the flux
            #ie, where q dot n is negative
            if dotprod < 0.0 and abs(dotprod) <= 1.0:
                q_i = abs(dotprod) * self.q_mag
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
    
    