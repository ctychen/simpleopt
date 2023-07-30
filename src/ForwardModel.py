import numpy as np

import Solid

def calculateIntegratedEnergy(q_mesh_all, trimeshSolid):
    faceAreas = trimeshSolid.area_faces
    mesh_q_dot_areas = q_mesh_all * faceAreas
    prods = np.where(mesh_q_dot_areas > 0, mesh_q_dot_areas, 0)
    mesh_energy = np.sum(np.abs(prods))
    return mesh_energy

def eich_profile(s_bar, q0, lambda_q, S, q_BG):
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


class ForwardModel_MeshHF: 

    def __init__(self, solidObj, q_mag, q_dir, hfMode): 
        self.solidObj = solidObj
        self.q_mag = q_mag #magnitude of applied q [W/m^2]
        self.q_dir = q_dir #direction of applied q, [x,y,z] [m]
        self.hfMode = hfMode #options: 'uniform', 'exponnorm', etc. 
        self.makeHFProfile(self.solidObj.trimeshSolid)
        # self.makeHFProfile(self.solidObj.trimeshSolid, self.q_dir)
        return
        
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
    
    