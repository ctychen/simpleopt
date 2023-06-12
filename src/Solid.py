import CADClass
import FreeCAD 
import Part
import Mesh
import MeshPart
import os
import numpy as np
import toolsClass
import math 

tools = toolsClass.tools()

class Box(CADClass.CAD): 
    
    # def __init__(self, stpfile):
    #     super(CADClass.CAD, self).__init__()
    #     self.STPfile = stpfile #want to be able to use as many HEAT cadclass fcns as possible, is this how to do it?
    #     return

    def __init__(self, stlfile, stpfile=""):
        super(CADClass.CAD, self).__init__()
        self.STLfile = stlfile
        self.STPfile = stpfile
        self.parts = self.loadSTEP()
        return

    #todo: maybe can replace some stuff wiht objFromPartnum from CADClass directly

    def createMesh(self, res=1000): 
        self.meshes = self.part2mesh(self.CADparts, res)
        # print("Meshed")
        return self.meshes

    def processModel(self, res=1000):
        #meshing the thing
        meshes = self.createMesh(res)
        if type(meshes) != list:
            meshes = [self.meshes]        
        #calc centers, normals, areas for meshed
        normcenterarea = self.normsCentersAreas(meshes)
        self.norms = normcenterarea[0] #norm[i] = [xi, yi, zi]
        self.centers = normcenterarea[1]
        self.areas = normcenterarea[2]
        print("Meshed, founds norms, centers, areas")
        return

    def rotateByAmount(self, thetax, thetay, thetaz, x, y, z):
        
        rot = FreeCAD.Placement(FreeCAD.Vector(x, y, z), FreeCAD.Rotation(thetax, thetay, thetaz))
        # print(len(FreeCAD.ActiveDocument.Objects))
        for obj in FreeCAD.ActiveDocument.Objects:
            if type(obj) == Part.Feature:
                print(f"Before Modifying Placement: {obj.Placement}")
                obj.Placement = rot.multiply(obj.Placement)
                obj.recompute()
                print(f"After Modifying Placement: {obj.Placement}")
        print("CAD Permutation Complete")
        return 

    def getCurrentRotation(self):
        return FreeCAD.ActiveDocument.Objects[0].Placement.Rotation #axis, angle

    def getCurrentRotationAngles(self):
        return self.getCurrentRotation().toEuler() #to access components, do (output)[0], etc. from 0-2
    
    def rotateTo(self, rotX, rotY, rotZ): #use this rotate fcn for now!
        rotMatrix = FreeCAD.Matrix() #unity matrix
        rotMatrix.rotateX(rotX)
        rotMatrix.rotateY(rotY)
        rotMatrix.rotateZ(rotZ)

        targetRotation = FreeCAD.Rotation(rotMatrix)

        for obj in FreeCAD.ActiveDocument.Objects:
            if type(obj) == Part.Feature:
                origPlacement = obj.Placement
                origRotation = origPlacement.Rotation
                relativeRotation = targetRotation.multiply(origRotation.inverted()) 

                newPlacement = origPlacement 
                newPlacement.Rotation = relativeRotation.multiply(origPlacement.Rotation)

                obj.Placement = newPlacement
                obj.recompute()
        print("Completed Rotation")
        return 



    def rotateModel(self, rotationMatrix4D):
        #FreeCAD.Placement can be (axis, angle) or as quaternion -> but for quaternion, is (x,y,z,w) --> axis, angle, or can be a Matrix4D. 
        #whatever is simpler

        rot = FreeCAD.Placement(rotationMatrix4D) #and then just use as same placement, ig -> can go between quaternion, Matrix4D in Rotation apparently?

        for obj in FreeCAD.ActiveDocument.Objects:
            if type(obj) == Part.Feature:
                print(f"Before Rotation: {obj.Placement}")
                obj.Placement = rot.multiply(obj.Placement)
                obj.recompute()
                print(f"After Rotation: {obj.Placement}")
        print("CAD Rotation Complete")
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
    

from scipy.spatial.transform import Rotation

class Box_Vector(CADClass.CAD): 
    
    def __init__(self, stlfile="", stpfile="", meshres=100):
        super(CADClass.CAD, self).__init__()
        self.STLfile = stlfile
        self.STPfile = stpfile
        self.parts = self.loadSTEP()

        self.allmeshes = self.part2meshStandard(self.CADparts) #self.part2mesh(self.CADparts, meshres) #this returns a list of meshes - replaced self.meshes but name confusion
        self.mesh = self.allmeshes[0] #this is a meshobject
        #print(f"Mesh: {self.mesh}")

        self.meshPoints = np.array(self.allmeshes[0].Points)

        self.verticesFromFacets = []

        for i in range(len(self.allmeshes[0].Facets)):
            facet_points = self.allmeshes[0].Facets[i].Points
            for point in facet_points:
                self.verticesFromFacets.append(list(point))
        
        # print(f"Vertices from facets: {self.verticesFromFacets}")
        # print(f"Length of vertices from points: {len(self.verticesFromFacets)}")

        self.initial_rotation = self.getCurrentRotationAngles()

        return
    

    def calculateRotationOnMesh(self, vertices, xAng, yAng, zAng): #is representing with angles even the way that makes the most sense???
        rotation = Rotation.from_euler('xyz', np.radians([xAng, yAng, zAng]), degrees=False)
        rotatedVertices = rotation.apply(vertices)

        return [rotation, rotatedVertices]
    
    def calculateRotationOnVector(self, vector, angles):
        #angles should be the target angle: rotate the vector to (angles), not by (angles)
        angles = np.radians(angles)
        r = Rotation.from_euler('xyz', angles)
        rotatedVector = r.apply(vector)        
        return rotatedVector
    
    def setVertices(self, newVertices):
        self.meshVertices = newVertices
        return 
    
    def findMeshCenter(self, mesh):
        center = mesh.BoundBox.Center
        return center
    
    
    def updateMesh(self, newVertices):
        if type(newVertices) != list:
            newVertices = newVertices.tolist()

        newMesh = Mesh.Mesh(newVertices)
        self.allmeshes[0] = newMesh
        return newMesh
    
    def makeMesh(self, vertices):
        if type(vertices) != list:
            vertices = vertices.tolist()

        newMesh = Mesh.Mesh(vertices)
        return newMesh
    
    def faceNormals(self, mesh):
        """
        returns normal vectors for single freecad mesh object in cartesian
        coordinates
        """
        #face normals
        normals = []
        for i, facet in enumerate(mesh.Facets):
            vec = np.zeros((3))
            for j in range(3):
                vec[j] = facet.Normal[j]
            normals.append(vec)
        return np.asarray(normals)

    def faceAreas(self, mesh):
        """
        returns face areas for mesh element
        """
        #face area
        areas = []
        for i, facet in enumerate(mesh.Facets):
            areas.append(facet.Area)
        return np.asarray(areas)


    def faceCenters(self, x, y, z):
        """
        returns centers of freecad mesh triangle in cartesian coordinates
        """
        #face centers
        centers = np.zeros((len(x), 3))
        centers[:,0] = np.sum(x,axis=1)/3.0
        centers[:,1] = np.sum(y,axis=1)/3.0
        centers[:,2] = np.sum(z,axis=1)/3.0
        return centers
    
    def calculateNormals(self, vertices):
        mesh = self.makeMesh(vertices)
        return self.faceNormals(mesh)
    
    def normsCentersAreas_VectorAny(self, vertices):
        mesh = self.makeMesh(vertices)
        norms = self.faceNormals(mesh)
        centers = [] #self.faceCenters(mesh)
        areas = [] #self.faceAreas(mesh) 
        return norms, centers, areas
    
    def getStandardMeshNorms(self): #use this for starting point
        standardMesh = self.part2meshStandard(self.CADparts)[0]
        normals = self.faceNormals(standardMesh)
        return normals
    
    def normsCentersAreas_Vector(self):
        return self.normsCentersAreas_VectorAny(self.meshVertices)
    
    def processModel(self): 
        #no need to re-mesh if we're doing vertices 
        # self.norms = self.calculateNorms(self.meshVertices)
        # self.centers = self.calculateCenters(self.meshVertices)
        # self.areas = self.calculateAreas(self.meshVertices)
        self.norms, self.centers, self.areas = self.normsCenterAreas_Vector()
        # print("Found norms, centers, areas for vectorized setup")
        return

    def getCurrentRotation(self):
        return FreeCAD.ActiveDocument.Objects[0].Placement.Rotation #axis, angle

    def getCurrentRotationAngles(self):
        return self.getCurrentRotation().toEuler() #to access components, do (output)[0], etc. from 0-2

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
    


class Box_Vector_Mesh(CADClass.CAD):

    def __init__(self, stlfile="", stpfile="", meshres=2):
        super(CADClass.CAD, self).__init__()
        self.STLfile = stlfile
        self.STPfile = stpfile
        self.parts = self.loadSTEP()

        self.allmeshes = self.part2mesh(self.CADparts, meshres)
        # self.allmeshes = self.part2meshStandard(self.CADparts)

        self.vertices = np.array(self.allmeshes[0].Points)

        self.faces = np.array(self.allmeshes[0].Facets)

        print(f"Vertices: {self.vertices}")

        print(f"Faces: {self.faces}")

        # self.allmeshes = self.part2meshStandard(self.CADparts) #self.part2mesh(self.CADparts, meshres) #this returns a list of meshes - replaced self.meshes but name confusion
        # self.mesh = self.allmeshes[0] #this is a meshobject
        #print(f"Mesh: {self.mesh}")

        # self.meshPoints = np.array(self.allmeshes[0].Points)

        # self.verticesFromFacets = []

        # for i in range(len(self.allmeshes[0].Facets)):
        #     facet_points = self.allmeshes[0].Facets[i].Points
        #     for point in facet_points:
        #         self.verticesFromFacets.append(list(point))
        
        # print(f"Vertices from facets: {self.verticesFromFacets}")
        # print(f"Length of vertices from points: {len(self.verticesFromFacets)}")

        # self.initial_rotation = self.getCurrentRotationAngles()

        return

    def makeMesh(self, vertices):
        if type(vertices) != list:
            vertices = vertices.tolist()

        newMesh = Mesh.Mesh(vertices)
        return newMesh
   
   #this one was prev attempts
    # def makePyramidFromCube(self, id=9):

    #     print(f"Starting pyramid attempt")

    #     vertices = []

    #     for i in range(len(self.allmeshes[0].Facets)):
    #         facet_points = self.allmeshes[0].Facets[i].Points
    #         for point in facet_points:
    #             vertices.append(list(point))      

    #     vertices = np.array(vertices)

    #     print(f"Vertices: {vertices}")  
    #     print(type(vertices))

    #     os.makedirs(f"pyramidtest{id}")
    #     meshUpdated = self.makeMesh(vertices)
    #     self.saveMeshSTL(meshUpdated, f"pyramidtest{id}/before_pyramid", "standard")

    #     faces = np.array([[facet.PointIndices[i] for i in range(3)] for facet in self.faces])
    #     print(f"Faces: {faces}")

    #     # base_facets = []
    #     # for facet in faces:
    #     #     if np.count_nonzero(np.isclose(vertices[facet][:, 2], -127.0)) == 3:
    #     #         base_facets.append(facet)
    #     # base_facets = np.array(base_facets)

    #     # Identify the bottom face of the cube
    #     bottom_face_indices = np.where(np.isclose(vertices[:, 2], 0))[0] #it's -127.0 right now bc that's the cube placement, but that should be standardized
    #     bottom_face = faces[np.isin(faces[:, 0], bottom_face_indices)]

    #     # Calculate the center of the bottom face
    #     center = np.mean(vertices[bottom_face[:, 1:]], axis=0)

    #     # Modify the vertices to form a pyramid
    #     for i in range(bottom_face.shape[0]):
    #         face_index = bottom_face[i, 0]
    #         vertex_indices = bottom_face[i, 1:]
    #         vertices[vertex_indices] = center

    #         # Update the face with the modified vertex indices
    #         faces[face_index, 1:] = vertex_indices
    #         meshUpdated = self.makeMesh(vertices)
    #         print(f"Making new mesh: {meshUpdated}")

    #         print(f"New vertices: {vertices}")
    #         print(f"New faces: {faces}")

    #         self.saveMeshSTL(meshUpdated, f"pyramidtest{id}/pyramid_test_{i}", "standard")

    #     # # Calculate the center of the bottom face
    #     # center = [sum(vertices[idx][i] for idx in bottom_face[0]) / 3 for i in range(3)]

    #     # # Create a mapping of old vertices to new vertices
    #     # vertex_map = {idx: idx for idx in range(len(vertices))}

    #     # # Modify the vertices to form a pyramid
    #     # count = 0
    #     # for face in bottom_face:
    #     #     for idx in face:
    #     #         if vertices[idx][2] != 0:
    #     #             vertices[idx] = center
    #     #             vertex_map[idx] = len(vertices)
            
    #     #     meshUpdated = self.makeMesh(vertices)
    #     #     print(f"Making new mesh: {meshUpdated}")

    #     #     print(f"New vertices: {vertices}")
    #     #     print(f"New faces: {faces}")

    #     #     self.saveMeshSTL(meshUpdated, f"pyramidtest{id}/pyramid_test_{count}", "standard")
    #     #     count += 1
            

    #     # Create a new list of modified faces with updated vertex indices
    #     # new_faces = []
    #     # for face in faces:
    #     #     new_face = [vertex_map[idx] for idx in face]
    #     #     new_faces.append(new_face)

    #     # faces = new_faces

    #     meshUpdated = self.makeMesh(vertices)
    #     print(f"Making new mesh: {meshUpdated}")

    #     print(f"New vertices: {vertices}")
    #     print(f"New faces: {faces}")

    #     self.saveMeshSTL(meshUpdated, f"pyramidtest{id}/pyramid_test_final", "standard")

    #     return vertices, faces, meshUpdated
    
    #for attempt at control-point-based process
    def compute_weights(self, vertices, control_points):
        distances = np.linalg.norm(vertices[:, np.newaxis] - control_points, axis=2)
        inv_distances = np.reciprocal(distances, where=distances != 0)
        weights = inv_distances / np.sum(inv_distances, axis=1, keepdims=True)
        return weights

    #use this one for now
    def pyramidFromCube(self, id=17):

        print(f"Starting pyramid attempt")

        vertices = []

        for i in range(len(self.allmeshes[0].Facets)):
            facet_points = self.allmeshes[0].Facets[i].Points
            for point in facet_points:
                vertices.append(list(point))      

        vertices = np.array(vertices)

        print(f"Original vertices: {vertices}")  
        print(f"Number of vertices: {len(vertices)}")

        os.makedirs(f"pyramidtest{id}")
        meshUpdated = self.makeMesh(vertices)
        self.saveMeshSTL(meshUpdated, f"pyramidtest{id}/before_pyramid", 2)

        faces = np.array([[facet.PointIndices[i] for i in range(3)] for facet in self.faces])
        print(f"Original faces: {faces}")    
        print(f"Number of faces: {faces.shape}")    

        control_points = np.array([
            [0, 10.0, 0],  #corner
            [10.0, 10.0, 0],   #corner
            [10.0, 0, 0],    #corner
            [0, 0, 0],   #corner
            #[0, 0, 10.0]      # Top Vertex
            [5.0, 5.0, 10], #top vertex

            #time to define some anchors along the edges to keep the base square - midpoints of base edges
            [5, 0, 0],
            [10, 5, 0],
            [5, 10, 0],
            [0, 5, 0]

        ])

        weights = self.compute_weights(vertices, control_points)

        # Deform the mesh by updating vertex positions based on control point weights
        deformed_vertices = np.zeros_like(vertices)
        count = 0
        for i in range(len(vertices)):
            deformed_vertices[i] = np.dot(weights[i], control_points)  
            meshUpdated = self.makeMesh(deformed_vertices)
            if count % 12 == 0: 
                self.saveMeshSTL(meshUpdated, f"pyramidtest{id}/pyramid_test_{i}", 2)
            count += 1

        meshUpdated = self.makeMesh(deformed_vertices)
        print(f"Making new mesh: {meshUpdated}")

        print(f"New vertices: {deformed_vertices}")
        #print(f"New faces: {faces}")
        print(f"Number of new vertices: {len(deformed_vertices)}")
        #print(f"Number of new faces: {len(faces)}")

        self.saveMeshSTL(meshUpdated, f"pyramidtest{id}/pyramid_test_final", 2)

        return deformed_vertices, faces
    
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
    
    