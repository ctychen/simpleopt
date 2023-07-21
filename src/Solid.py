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

### Mesh Properties ###

def calculateFaceNormals(vertices, faces):
    """
    calculate normal vectors for each face in a trimesh.

    Parameters
    ----------
    vertices : (n, 3) array
        Vertex coordinate in the space.
    faces : (m, 3) array
        Indices of vertices that make up each face.
    
    Returns
    -------
    face_normals : (m, 3) array
        Normal vector for each face.
    """

    # #get vectors of the two edges for each face
    # vec1 = vertices[faces[:, 1]] - vertices[faces[:, 0]]
    # vec2 = vertices[faces[:, 2]] - vertices[faces[:, 0]]
    
    # #calculate normal vectors using cross product
    # faceNormals = np.cross(vec1, vec2)
    
    # # Normalize each normal vector
    # norms = np.linalg.norm(faceNormals, axis=1)
    # faceNormals = faceNormals / norms[:, np.newaxis]

    #get vectors of the two edges for each face
    vecs = vertices[faces[:, 1:]] - vertices[faces[:, 0, np.newaxis]]
    
    #calculate normal vectors using cross product
    faceNormals = np.cross(vecs[:, 0], vecs[:, 1])
    
    # Normalize each normal vector
    norms = np.linalg.norm(faceNormals, axis=1, keepdims=True)
    faceNormals /= norms
    
    return faceNormals


def calculateVertexDefects(vertices, faces):
    """
    Compute the vertex defects for each vertex in a mesh - vertex defects 2*pi - sum of angles around each vertex

    Parameters
    ----------
    vertices : (n, 3) array
        Vertex coordinate in the space.
    faces : (m, 3) array
        Indices of vertices that make up each face.
    
    Returns
    -------
    vertex_defects : (n,) array
        Vertex defect for each vertex.
    """

    vertex_defects = 2 * np.pi * np.ones(vertices.shape[0])

    for face in faces:
        #vectors to previous and next vertices
        prev_vectors = vertices[face] - vertices[face[[2, 0, 1]]]
        next_vectors = vertices[face[[1, 2, 0]]] - vertices[face]

        #normalize vectors
        prev_vectors /= np.linalg.norm(prev_vectors, axis=-1, keepdims=True)
        next_vectors /= np.linalg.norm(next_vectors, axis=-1, keepdims=True)

        #angle at each vertex is the angle between the two vectors
        cosines = np.sum(prev_vectors * next_vectors, axis=-1)
        angles = np.arccos(np.clip(cosines, -1, 1))

        #subtract these angles from the vertices' defects
        np.subtract.at(vertex_defects, face, angles)

    return #vertex_defects


def calculateFaceCenters(vertices, faces):
    """
    Compute the center of each face.

    Args:
    vertices: A (n, 3) array representing the coordinates of each vertex.
    faces: A (m, 3) array representing the faces, where each row contains the indices of three vertices.

    Returns:
    A (m, 3) array representing the center of each face.
    """
    #index into the vertices array to get the vertices of each face.
    face_vertices = vertices[faces]
    
    #compute the center of each face by taking the mean across the second dimension.
    face_centers = np.mean(face_vertices, axis=1)

    return face_centers


def calculateFaceAreas(vertices, faces):
    """
    Compute the area of each face.

    Args:
    vertices: A (n, 3) array representing the coordinates of each vertex.
    faces: A (m, 3) array representing the faces, where each row contains the indices of three vertices.

    Returns:
    A (m,) array representing the area of each face.
    """
    #get vectors of the two edges for each face
    vec1 = vertices[faces[:, 1]] - vertices[faces[:, 0]]
    vec2 = vertices[faces[:, 2]] - vertices[faces[:, 0]]
    
    #calculate normal vectors using cross product
    face_normals = np.cross(vec1, vec2)
    
    #area of each face is half the norm of the cross product of the two edges
    face_areas = 0.5 * np.linalg.norm(face_normals, axis=1)
    
    return face_areas


class MeshSolid(CADClass.CAD):

    def __init__(self, stlfile="", stpfile="", meshres=2.0):
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
        save vertices, faces, and face adjacency to class variables - should take these in for calculation for normals, centers, areas, angle diffs, etc. 
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

        self.meshVertices = np.array(trimeshSolid.vertices) 
        self.meshFaces = np.array(trimeshSolid.faces)
        self.meshFaceAdjacency = np.array(trimeshSolid.face_adjacency)
        # self.meshNormals = Solid.calculateFaceNormals(self.meshVertices, self.meshFaces)

        return trimeshSolid    


    ### Old functions ###

    # def normalsCentersAreas_Trimesh(self, trimeshSolid):

    #     """
    #     Return lists of normal vectors, centers, areas for each mesh face
    #     """

    #     #normal vectors of each face
    #     faceNormals = trimeshSolid.face_normals

    #     #center of each face
    #     faceCenters = trimeshSolid.triangles_center

    #     #area of each face
    #     faceAreas = trimeshSolid.area_faces

    #     return faceNormals, faceCenters, faceAreas


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

    
    