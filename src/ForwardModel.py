import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

class ForwardModel_MeshHF: 

    def __init__(self, solidObj, q_mag, q_dir): 
        self.solidObj = solidObj
        self.q_mag = q_mag #magnitude of applied q [W/m^2]
        self.q_dir = q_dir #direction of applied q, [x,y,z] [m]
        return
    
    
    #overall: can't have HF on backface - and doesn't make sense to have negative HF. 
    #so should be setting backfaces (negative HF) to 0
    #this might also eliminate some weirdness with back face corners/edges moving
    #calc: -q*(b dot n)
    #then take any values where -(-q*(b dot n)) > 0, ie. HF's aren't physical
    #and set those to 0

    def calculateAllHF_AllVals(self, trimeshSolid):
        """
        Calculate HF on every mesh element, with no exclusions, and returning them all in a list
        """
        
        normals = trimeshSolid.face_normals
        q_mesh_all = -1 * (np.dot(normals, self.q_dir)) * self.q_mag
        return q_mesh_all
    

    def calculateAllHF(self, trimeshSolid):
        """
        Calculate HF on every mesh element, but set negative values to 0 (remove nonphysical element values)
        """
        q_mesh_all = self.calculateAllHF_AllVals(trimeshSolid)

        return q_mesh_all
    
    
    def calculateHFDistribution(self, trimeshSolid):
        """
        Calculate mean, variance, standard deviation on all HF over mesh elements
        """
        
        normals = trimeshSolid.face_normals
        q_mesh_all = -1 * (np.dot(normals, self.q_dir)) * self.q_mag

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

        normals = trimeshSolid.face_normals
        q_mesh_all = -1 * (np.dot(normals, self.q_dir)) * self.q_mag

        dotprods = filter(lambda x: x > 0, q_mesh_all)
        q_mesh_sum = sum(abs(x) for x in dotprods)

        return q_mesh_sum
    
    def calculateIntegratedEnergy(self, trimeshSolid):
        normals = trimeshSolid.face_normals
        faceAreas = np.array(trimeshSolid.area_faces)
        q_vals = np.array(-1 * (np.dot(normals, self.q_dir)) * self.q_mag)
        mesh_q_dot_areas = q_vals * faceAreas
        prods = filter(lambda x: x > 0, mesh_q_dot_areas)
        mesh_energy = sum(abs(x) for x in prods)
        
        return mesh_energy
    
    def calculateMaxHF(self, trimeshSolid):
        """
        Find single highest heat flux from all mesh element heat fluxes
        """
        q_mesh_all = self.calculateAllHF(trimeshSolid)
        return np.max(q_mesh_all)
    
    
    