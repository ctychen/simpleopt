import numpy as np

class FowardModel_Box: 
    def __init__(self, solid, q_mag, q_dir, g_x): #some arbitrary starting values todo
        self.solid = solid
        self.q_mag = q_mag #magnitude of applied q [W/m^2]
        self.q_dir = q_dir #direction of applied q, [x,y,z] [m]
        self.g_x = g_x #objective function, can swap out
        return
    
    def processCADModel(self):
        #meshing the thing
        self.meshes = self.solid.createMesh()
        #calc centers, normals, areas for meshed
        normcenterarea = self.solid.normsCentersAreas(self.meshes)
        self.norms = normcenterarea[0] #norm[i] = [xi, yi, zi]
        self.centers = normcenterarea[1]
        self.areas = normcenterarea[2]
        return
    
    def calcQMesh(self):
        #for every mesh element, calculate q and store in list/arr with same length as mesh elem list
        q_mesh_all = []
        for i in range(len(self.meshes)): 
            #calculate q
            #defn for q: (qdir dot normal)*(mag of q) at center for each element
            #how to handle q? should be uniform so direction and mag constant
            #add q value to list
            norm = self.norms[i]
            q_i = np.dot(self.q_dir, norm) * self.q_mag
            q_mesh_all.append(q_i)

        return q_mesh_all
    
    def calcObjective(self, q_mesh_all): 
        #takes in objective function
        return self.g_x(q_mesh_all)
    