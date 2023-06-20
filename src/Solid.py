import CADClass
import FreeCAD 
import Part
import Mesh
import MeshPart
import os
import numpy as np
import toolsClass
import math 

import trimesh
# Check available backends
print(trimesh.interfaces.blender.exists)
print(trimesh.interfaces.scad.exists)

tools = toolsClass.tools()

class Box_Vector_Mesh(CADClass.CAD):

    def __init__(self, stlfile="", stpfile="", meshres=1.0):
        super(CADClass.CAD, self).__init__()
        self.STLfile = stlfile
        self.STPfile = stpfile
        self.parts = self.loadSTEP()

        self.meshres = meshres

        self.allmeshes = self.part2mesh(self.CADparts, meshres)
        # self.allmeshes = self.part2meshStandard(self.CADparts)

        self.vertices = np.array(self.allmeshes[0].Points)
        self.faces = np.array(self.allmeshes[0].Facets)

        return

    def makeMesh(self, vertices):
        """
        create freecad mesh from array or list of vertices
        """
        if type(vertices) != list:
            vertices = vertices.tolist()

        newMesh = Mesh.Mesh(vertices)
        return newMesh
    
    def gradMoveVertex_TriMesh(self, vertex, tri_mesh, step_size = 0.1):
        """
        gradient descent implementation for trimesh mesh
        gradually move a vertex in the direction such that it approaches the closest point on the target mesh
        """
        # closest_point, distance, triangle_id = mesh.nearest.on_surface([point])
        closest_point, distance, triangle_id = tri_mesh.nearest.on_surface([vertex])
        # print(f"Closest point on the mesh: {closest_point[0]}")
        # print(f"Distance from point to mesh: {distance[0]}")
        direction_vector = closest_point[0] - vertex
        #if norms != 0.0: #so we don't have division by 0 issues
        #    direction_vector /= np.linalg.norm(direction_vector)
        #move vertex towards closest point using direction vector
        vertex += step_size * direction_vector
        return vertex
    

    def gradientDescentHF(self, tri_mesh, objectiveFunction, delta):
        """
        gradient descent implementation for heat flux minimization
        takes in trimesh object, applies objective function to it
        moves all mesh vertices by a small amount and finds the gradient when doing so 
        """
        gradient = np.zeros_like(tri_mesh.Vertices)
        #for each vertex
        for i in range(len(tri_mesh.Vertices)):
            #for each dimension
            for j in range(3):
                #move vertex a bit, see if the objective function decreases
                tri_mesh.Vertices[i, j] += delta
                obj_afterMoving = objectiveFunction(tri_mesh)
                #move vertex back to original, calc objective function then for comparison
                tri_mesh.Vertices[i, j] -= delta
                obj_beforeMoving = objectiveFunction(tri_mesh)
                gradient[i, j] = (obj_afterMoving - obj_beforeMoving) / (2 * delta)
        return gradient


    def meshHFOpt(self, id='000'):
        """
        process to optimize surface heat flux using gradient descent
        """

        print(f"Starting HF optimization attempt")

        vertices = []

        for i in range(len(self.allmeshes[0].Facets)):
            facet_points = self.allmeshes[0].Facets[i].Points
            for point in facet_points:
                vertices.append(list(point))      

        vertices = np.array(vertices) #cube vertices now defined

        print(f"Number of vertices: {len(vertices)}")

        os.makedirs(f"pyramidtest{id}")
        meshUpdated = self.makeMesh(vertices)
        self.saveMeshSTL(meshUpdated, f"pyramidtest{id}/before_pyramid", self.meshres)

        faces = np.array([[facet.PointIndices[i] for i in range(3)] for facet in self.faces])
        # print(f"Original faces: {faces}")    
        print(f"Number of faces: {faces.shape}")   

        cubeVertices = np.array([[v.x, v.y, v.z] for v in self.allmeshes[0].Points])
        cubeFaces = np.array([[f.PointIndices[0], f.PointIndices[1], f.PointIndices[2]] for f in self.allmeshes[0].Facets])
        count = 0

        cubeTriMesh = trimesh.Trimesh(
            vertices=cubeVertices, faces=cubeFaces
        )

        cubeTriMesh.export(f"pyramidtest{id}/basecube.stl")
        
        def hfObjectiveFunction(cubeTrimesh):
            #TODO: techically already exists elsewhere, but gotta trimeshify
            #Calculate normal vector for every mesh facet
            #Calculate (dot product of normal vector, direction of HF)*(magnitude of HF)
            #Store all these in a list, then find the max
            meshAllHF = []
            return meshMaxHF
        

        while hfObjectiveFunction(cubeTriMesh) > 0:
            #calculated the objective function value and it means HF is still not @ min
            #move each vertex in direction to minimize HF - how to determine this/pick direction?  

            #calc the gradient
            gradient = self.gradientDescentHF()

            #move the vertices a bit based on the gradient (not sure if you can do this without looping)
            cubeTriMesh.Vertices -= 0.01 * gradient
        
        mesh2 = self.makeMesh(vertices)

        self.saveMeshSTL(mesh2, f"pyramidtest{id}/pyramid_test_final_{id}", self.meshres)

        return 
    

    def saveMeshSTL(self, mesh, label, resolution, path='./'):
        """
        Writes a mesh object to STL file named by part number.
        If mesh is a list of mesh objects, then write a separate file for
        each mesh object in the list.  Clobbers if overWriteMask is True
        """
        #Check if this is a single file or list and make it a list
        if type(mesh) != list:
            mesh = [mesh]
        if type(label)!= np.ndarray:
            if type(label) != list:
                label = [label]
        if type(resolution) != list:
            resolution=[resolution]*len(mesh)

        #Recursively make dirs for STLs
        print("making STL directory")
        tools.makeDir(path, clobberFlag=False, mode=0o774, UID=-1, GID=-1)

        for i in range(len(mesh)):
            # ___ (3 underdashes) is the str we use to separate mesh name from resolution
            # this MATTERS when we read back in a mesh (see self.loadROIMesh and self.loadIntersectMesh)

            #standard meshing algorithm
            stdList = ['standard', 'Standard', 'STANDARD']
            if resolution[i] in stdList:
                filename = path + label[i] + "___"+resolution[i]+".stl"
            #mefisto meshing algorithm
            else:
                filename = path + label[i] + "___{:.2f}mm.stl".format(float(resolution[i]))
            if os.path.exists(filename): # and self.overWriteMask == False:
                print("Not clobbering mesh file...")
            else:
                print("Writing mesh file: " + filename)
                mesh[i].write(filename)
                os.chmod(filename, 0o774)
                # os.chown(filename, self.UID, self.GID)

        print("\nWrote meshes to files")
        return
    
    