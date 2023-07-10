import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

class ForwardModel_MeshHF: 

    def __init__(self, solidObj, q_mag, q_dir, hfMode): 
        self.solidObj = solidObj
        self.q_mag = q_mag #magnitude of applied q [W/m^2]
        self.q_dir = q_dir #direction of applied q, [x,y,z] [m]
        self.hfMode = hfMode #options: 'uniform', 'exponnorm', etc. 
        self.makeHFProfile(self.solidObj.trimeshSolid)
        # self.makeHFProfile(self.solidObj.trimeshSolid, self.q_dir)
        return
    
    def eich_profile(self, s_bar, q0, lambda_q, S, q_BG):
        """
        Fits an Eich profile to heat flux data, using upstream-mapped
        distance to remove the explicit flux-expansion factor

        See Eich Nuclear Fusion 2013. Using s_bar as the OMP-mapped distance
        means that you should set 'f_x' = 1

        You can pass the parameters to the profile either all as dimensional
        Quantity objects, or all as dimensionless raw arrays/floats

        s_bar: OMP-mapped radial distance [mm]
        q0: peak heat flux [MW/m^2]
        lambda_q: heat flux decay width [mm]
        S: spreading factor [mm]
        q_BG: background heat-flux [MW/m^2]
        """

        from scipy.special import erfc

        return (
            q0/ 2.0
            * np.exp((S / (2.0 * lambda_q)) ** 2 - s_bar / lambda_q)
            * erfc(S / (2.0 * lambda_q) - s_bar / S)
            + q_BG
        )
    
    
    def makeHFProfile(self, trimeshSolid): #, directionVectors):
        """
        make profile with direction and magnitude of HF at each face center
        this fcn is bc we currently have a uniform direction of HF everywhere so just 1 vector ok, but if not, 
        use version above
        """
        meshHFProfile = []
        directionVectors = self.q_dir

        if self.hfMode == 'uniform':
            meshHFProfile = self.makeUniformHFProfile(trimeshSolid, directionVectors)
        elif self.hfMode == "uniform_multiple":
            meshHFProfile = self.makeMultipleHFProfile(trimeshSolid, directionVectors)
        elif self.hfMode == 'exponnorm':
            #but also we shouldn't be doing this for all faces, only top faces? 
            all_HF_magnitudes = self.calculateHFProfileMagnitudes(trimeshSolid)
            # for i in range(len(all_HF_magnitudes)):
            #     magnitude = all_HF_magnitudes[i]
            #     meshHFProfile.append([magnitude, directionVectors])

            for magnitude in all_HF_magnitudes:
                meshHFProfile.append([magnitude, directionVectors])

        self.meshHFProfile = meshHFProfile
        return meshHFProfile
    
    
    def makeUniformHFProfile(self, trimeshSolid, directionVectors):
        """ 
        return list of HF magnitudes and direction vectors for each face center
        """
        return [[self.q_mag, directionVectors]] * len(trimeshSolid.face_normals)
    
    def makeMultipleHFProfile(self, trimeshSolid, directionVectors):
        """
        for now, make profile with 2 different HF magnitudes and direction vectors - could extend to more eventually
        """
        # directions_magnitudes_list = [[mag, dir_vec] for mag, dir_vec in zip(self.q_mag, directionVectors)]
        # return directions_magnitudes_list * len(trimeshSolid.face_normals)
        return [[self.q_mag[0], directionVectors[0]], [self.q_mag[1], directionVectors[1]]] * len(trimeshSolid.face_normals)

    def calculateHFProfileMagnitudes(self, trimeshSolid):
        """
        calc magnitudes of nonuniform heat flux on surface as function of x-coordinate of face center
        eventually could use proper eich profile but for now using exponentially modified gaussian
        """
        from scipy.special import erfc
        from scipy.stats import exponnorm

        centers = trimeshSolid.triangles_center

        x_centers = centers[:,0] #isolate x-values of the centers

        # Define parameters
        K = 1.0  # This defines the "shape" of the curve (larger K -> more skew)
        mu = 5.0  # This is the mean of the normal part of the distribution
        sigma = 2.0 #1.0 #2.5 #1.0  # This is the standard deviation of the normal part
        q_mag_max = self.q_mag #for now use qmag as maximum for qmag distribution

        # Calculate the PDF at these x values
        q_mag_all_centers = exponnorm.pdf(x_centers, K, loc=mu, scale=sigma) 
        q_mag_all_centers = q_mag_all_centers * q_mag_max / np.max(q_mag_all_centers)

        # self.q_mag_all_centers = q_mag_all_centers
        return q_mag_all_centers
    
    
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
        calculate HF on every mesh element, but this should use nonuniform HF on mesh element centers
        this could also use uniform HF; just need meshHFProfile to be a list of the same HF value for each element
        """

        q_mesh_all = []

        if self.hfMode == 'uniform':
            q_mesh_all = self.calculateAllHF_AllVals(trimeshSolid)

        elif self.hfMode == 'exponnorm':
            normals = trimeshSolid.face_normals
            # newHFProfile = self.calculateHFProfileMagnitudes(trimeshSolid)
            newHFProfile = self.meshHFProfile
            for i in range(len(newHFProfile)):
                magnitude = newHFProfile[i][0]
                direction = newHFProfile[i][1]
                q_mesh_all.append(-1 * (np.dot(normals[i], direction)) * magnitude)
            # for i in range(len(self.meshHFProfile)):
            #     # magnitude = self.meshHFProfile[i][0]
            #     direction = self.meshHFProfile[i][1]
            #     q_mesh_all.append(-1 * (np.dot(normals[i], direction)) * magnitude)
            q_mesh_all = np.array(q_mesh_all)

        elif self.hfMode == "uniform_multiple":
            normals = trimeshSolid.face_normals
            # q_mesh_all = [
            #     sum(max(-1 * self.q_mag[j] * np.dot(n, self.q_dir[j]), 0) for j in range(len(self.q_dir))) 
            #     for n in normals
            # ]
            # q_mesh_all = np.array(q_mesh_all)
            for i in range(len(normals)):
                n = normals[i]
                q_i = 0
                for j in range(len(self.q_dir)):
                    q_to_add = (-1 * self.q_mag[j] * (np.dot(n, self.q_dir[j]))) if (-1 * np.dot(n, self.q_dir[j])) > 0 else 0
                    q_i += q_to_add
                q_mesh_all.append(q_i)
            q_mesh_all = np.array(q_mesh_all)

        return q_mesh_all
        

    def calculateHFMeshSum(self, q_mesh_all):
        """
        Calculate sum of heat flux from all mesh elements
        """

        q_mesh_vals = np.where(q_mesh_all > 0, q_mesh_all, 0)
        q_mesh_sum = np.sum(np.abs(q_mesh_vals))

        return q_mesh_sum
    
    
    def calculateIntegratedEnergy(self, q_mesh_all, trimeshSolid):
        faceAreas = trimeshSolid.area_faces
        mesh_q_dot_areas = q_mesh_all * faceAreas
        prods = np.where(mesh_q_dot_areas > 0, mesh_q_dot_areas, 0)
        mesh_energy = np.sum(np.abs(prods))
        return mesh_energy
    
    
    def calculateMaxHF(self, q_mesh_all):
        return np.max(q_mesh_all)
    

    def filteredCalculateMaxHF(self, q_mesh_all, unconstrainedFaces = []):
        """
        filtering out unconstrained faces from max HF calculation - originally was added bc of HF with incident angle
        if no unconstrained faces defined, then we just return max HF same as before
        """

        # if unconstrainedFaces: #if unconstrainedFaces is not empty
        #     unconstrainedFaces = list(unconstrainedFaces) 
        #     mask = np.ones(q_mesh_all.shape, dtype=bool)
        #     mask[unconstrainedFaces] = False
        #     q_mesh_all[mask] = 0

        return np.max(q_mesh_all)
    
    
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


    # def calculateHFMeshSum(self, trimeshSolid):
    #     """
    #     Calculate sum of heat flux from all mesh elements
    #     """

    #     normals = trimeshSolid.face_normals
    #     q_mesh_all = []
    #     for i in range(len(self.meshHFProfile)):
    #         magnitude = self.meshHFProfile[i][0]
    #         direction = self.meshHFProfile[i][1]
    #         q_mesh_all.append(-1 * (np.dot(normals[i], direction)) * magnitude)
    #     q_mesh_all = np.array(q_mesh_all)

    #     q_mesh_vals = np.where(q_mesh_all > 0, q_mesh_all, 0)
    #     q_mesh_sum = np.sum(np.abs(q_mesh_vals))

    #     return q_mesh_sum

    # def calculateIntegratedEnergy(self, trimeshSolid):
    #     normals = trimeshSolid.face_normals
    #     faceAreas = trimeshSolid.area_faces
    #     q_mesh_all = []
    #     for i in range(len(self.meshHFProfile)):
    #         magnitude = self.meshHFProfile[i][0]
    #         direction = self.meshHFProfile[i][1]
    #         q_mesh_all.append(-1 * (np.dot(normals[i], direction)) * magnitude)
    #     q_mesh_all = np.array(q_mesh_all)

    #     mesh_q_dot_areas = q_mesh_all * faceAreas

    #     prods = np.where(mesh_q_dot_areas > 0, mesh_q_dot_areas, 0)
    #     mesh_energy = np.sum(np.abs(prods))

    #     return mesh_energy
    
    # def calculateMaxHF(self, trimeshSolid):
    #     """
    #     Find single highest heat flux from all mesh element heat fluxes
    #     """
    #     q_mesh_all = self.calculateAllHF(trimeshSolid)
    #     return np.max(q_mesh_all)


    
    
    