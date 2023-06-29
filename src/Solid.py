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

class MeshSolid(CADClass.CAD):

    def __init__(self, stlfile="", stpfile="", meshres=2.5):
        super(CADClass.CAD, self).__init__()
        self.STLfile = stlfile
        self.STPfile = stpfile
        self.parts = self.loadSTEP()

        self.meshres = meshres

        self.allmeshes = self.part2mesh(self.CADparts, meshres)
        # self.allmeshes = self.part2meshStandard(self.CADparts)

        self.vertices = np.array(self.allmeshes[0].Points)
        self.faces = np.array(self.allmeshes[0].Facets)

        self.trimeshSolid = self.processSolid()

        return

    
    def makeFreecadMesh(self, vertices):
        """
        create freecad mesh from array or list of vertices
        """
        if type(vertices) != list:
            vertices = vertices.tolist()

        newMesh = Mesh.Mesh(vertices)
        return newMesh

    def processSolid(self):
        """
        take solid, mesh it, and convert it into a trimesh mesh object
        """

        vertices = []

        for i in range(len(self.allmeshes[0].Facets)):
            facet_points = self.allmeshes[0].Facets[i].Points
            for point in facet_points:
                vertices.append(list(point))      

        vertices = np.array(vertices) #cube vertices now defined

        print(f"Number of vertices: {len(vertices)}")

        faces = np.array([[facet.PointIndices[i] for i in range(3)] for facet in self.faces])
        # print(f"Original faces: {faces}")    
        print(f"Number of faces: {faces.shape}")   

        solidVertices = np.array([[v.x, v.y, v.z] for v in self.allmeshes[0].Points])
        solidFaces = np.array([[f.PointIndices[0], f.PointIndices[1], f.PointIndices[2]] for f in self.allmeshes[0].Facets])

        trimeshSolid = trimesh.Trimesh(
            vertices=solidVertices, faces=solidFaces
        )     

        return trimeshSolid
    

    def normalsCentersAreas_Trimesh(self, trimeshSolid):

        """
        Return lists of normal vectors, centers, areas for each mesh face
        """

        #normal vectors of each face
        faceNormals = trimeshSolid.face_normals

        #center of each face
        faceCenters = trimeshSolid.triangles_center

        #area of each face
        faceAreas = trimeshSolid.area_faces

        return faceNormals, faceCenters, faceAreas


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

    
    