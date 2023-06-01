import numpy as np

class ForwardModel_Box: 
    def __init__(self, g_x, solid, q_mag, q_dir): #some arbitrary starting values todo
        self.solid = solid
        self.q_mag = q_mag #magnitude of applied q [W/m^2]
        self.q_dir = q_dir #direction of applied q, [x,y,z] [m]
        self.g_x = g_x #objective function, can swap out
        return
    
    # todo: maybe this shouldn't stay here. and instead be in solid.py
    def processCADModel(self):
        #meshing the thing
        meshes = self.solid.createMesh(res=2000)
        if type(meshes) != list:
            meshes = [self.meshes]        
        #calc centers, normals, areas for meshed
        normcenterarea = self.solid.normsCentersAreas(meshes)
        self.norms = normcenterarea[0] #norm[i] = [xi, yi, zi]
        self.centers = normcenterarea[1]
        self.areas = normcenterarea[2]
        # print("Model processed")
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
    
    def calcObjective(self, q_mesh_all): 
        #takes in objective function
        return self.g_x(q_mesh_all)
    