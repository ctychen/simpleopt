import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

class ForwardModel_MeshHF: 

    def __init__(self, solidObj, q_mag, q_dir): 
        self.solidObj = solidObj
        self.q_mag = q_mag #magnitude of applied q [W/m^2]
        self.q_dir = q_dir #direction of applied q, [x,y,z] [m]

        return
    
    
    def makeHFProfile(self, trimeshSolid, directionVector):
        """
        make profile with direction and magnitude of HF at each face center
        this fcn is bc we currently have a uniform direction of HF everywhere so just 1 vector ok, but if not, 
        use version above
        """
        all_center_HFs = []
        #but also we shouldn't be doing this for all faces, only top faces? 
        all_HF_magnitudes = self.calculateHFProfileMagnitudes(trimeshSolid)
        # for i in range(len(all_HF_magnitudes)):
        #     magnitude = all_HF_magnitudes[i]
        #     all_center_HFs.append([magnitude, directionVector])

        for magnitude in all_HF_magnitudes:
            all_center_HFs.append([magnitude, directionVector])

        # print(f"Made HF profile: {all_center_HFs}")
        self.all_center_HFs = all_center_HFs
        return all_center_HFs
    
    
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
    

    # def calculateAllHF(self, trimeshSolid):
    #     """
    #     Calculate HF on every mesh element
    #     """
    #     q_mesh_all = self.calculateAllHF_AllVals(trimeshSolid)

    #     return q_mesh_all
    
    def calculateAllHF(self, trimeshSolid):
        """
        calculate HF on every mesh element, but this should use nonuniform HF on mesh element centers
        """
        normals = trimeshSolid.face_normals
        #self.all_center_HFs is a list of lists, where each sublist is [magnitude, direction]
        #so we want to get the magnitudes and multiply that by dot product of direction vector with normal vector like before
        #self.all_center_HFs[i][0] is the magnitude of the ith element, and self.all_center_HFs[i][1] is the direction vector
        #so we want to get the dot product of self.all_center_HFs[i][1] with normals
        #and then multiply that by self.all_center_HFs[i][0]
        
        q_mesh_all = []
        for i in range(len(self.all_center_HFs)):
            magnitude = self.all_center_HFs[i][0]
            direction = self.all_center_HFs[i][1]
            q_mesh_all.append(-1 * (np.dot(normals[i], direction)) * magnitude)
        q_mesh_all = np.array(q_mesh_all)
        # print(f"q_mesh_all: {q_mesh_all}")
        return q_mesh_all
    

    def calculateHFProfileMagnitudes(self, trimeshSolid):
        """
        calc magnitudes of nonuniform heat flux on surface as function of x-coordinate of face center
        eventually could use proper eich profile but for now using exponentially modified gaussian
        """
        from scipy.special import erfc
        centers = trimeshSolid.triangles_center
        
        # #want to get x-vals of centers
        # mesh_center_yvals = centers[:, 1]
        # unconstrainedFaces = list(set(np.where(mesh_center_yvals == 10.0)[0])) #this should isolate the y-values
        # centers = centers[np.array(unconstrainedFaces)] #keep only indices of faces on top surface

        x_centers = centers[:,0] #isolate x-values of the centers

        #for this overall we're calculating HF for every mesh element since this way don't need to track more indices for HF calculation later 

        mean = 5.0
        sigma = 1.0
        l = 1.0
        q_mag_max = self.q_mag #for now use qmag as maximum for qmag distribution

        # q_all = 0.5 * np.exp(rho_0**2 - rho) * erfc(rho_0 - rho/(2*rho_0))
        q_mag_all_centers = (q_mag_max) * (l / 2) * np.exp((l/2)*(2*mean + l*sigma**2 - 2*x_centers)) * erfc((mean + l*sigma**2 - x_centers)/(np.sqrt(2) * sigma))
        # print(f"q magnitude on centers: {q_mag_all_centers}")
        # input()
        self.q_mag_all_centers = q_mag_all_centers
        return q_mag_all_centers
    
    # def calculateHFMeshSum(self, trimeshSolid):
    #     """
    #     Calculate sum of heat flux from all mesh elements
    #     """

    #     normals = trimeshSolid.face_normals
    #     q_mesh_all = -1 * (np.dot(normals, self.q_dir)) * self.q_mag

    #     q_mesh_vals = np.where(q_mesh_all > 0, q_mesh_all, 0)
    #     q_mesh_sum = np.sum(np.abs(q_mesh_vals))

    #     return q_mesh_sum

    def calculateHFMeshSum(self, trimeshSolid):
        """
        Calculate sum of heat flux from all mesh elements
        """

        normals = trimeshSolid.face_normals
        q_mesh_all = []
        for i in range(len(self.all_center_HFs)):
            magnitude = self.all_center_HFs[i][0]
            direction = self.all_center_HFs[i][1]
            q_mesh_all.append(-1 * (np.dot(normals[i], direction)) * magnitude)
        q_mesh_all = np.array(q_mesh_all)

        q_mesh_vals = np.where(q_mesh_all > 0, q_mesh_all, 0)
        q_mesh_sum = np.sum(np.abs(q_mesh_vals))

        return q_mesh_sum
    

    # def calculateIntegratedEnergy(self, trimeshSolid):
    #     normals = trimeshSolid.face_normals
    #     faceAreas = trimeshSolid.area_faces
    #     q_vals = -1 * (np.dot(normals, self.q_dir)) * self.q_mag
    #     mesh_q_dot_areas = q_vals * faceAreas

    #     prods = np.where(mesh_q_dot_areas > 0, mesh_q_dot_areas, 0)
    #     mesh_energy = np.sum(np.abs(prods))

    #     return mesh_energy

    def calculateIntegratedEnergy(self, trimeshSolid):
        normals = trimeshSolid.face_normals
        faceAreas = trimeshSolid.area_faces
        q_mesh_all = []
        for i in range(len(self.all_center_HFs)):
            magnitude = self.all_center_HFs[i][0]
            direction = self.all_center_HFs[i][1]
            q_mesh_all.append(-1 * (np.dot(normals[i], direction)) * magnitude)
        q_mesh_all = np.array(q_mesh_all)

        mesh_q_dot_areas = q_mesh_all * faceAreas

        prods = np.where(mesh_q_dot_areas > 0, mesh_q_dot_areas, 0)
        mesh_energy = np.sum(np.abs(prods))

        return mesh_energy
    
    
    def calculateMaxHF(self, trimeshSolid):
        """
        Find single highest heat flux from all mesh element heat fluxes
        """
        q_mesh_all = self.calculateAllHF(trimeshSolid)
        return np.max(q_mesh_all)


    
    # def makeHFProfile(self, trimeshSolid, all_directions):
    #     """
    #     make profile with direction and magnitude of HF at each face center
    #     """
    #     all_center_HFs = []
    #     all_HF_magnitudes = self.calculateHFProfileMagnitudes(trimeshSolid)
    #     for i in range(len(all_HF_magnitudes)):
    #         magnitude = all_HF_magnitudes[i]
    #         direction = all_directions[i]
    #         all_center_HFs.append([magnitude, direction])
    #     print(f"Made HF profile: {all_center_HFs}")
    #     return all_center_HFs
    
    
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
    
    
    
    