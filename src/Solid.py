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
from scipy.optimize import minimize

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

        # print(f"Vertices: {self.vertices}")

        # print(f"Faces: {self.faces}")

        return

    def makeMesh(self, vertices):
        if type(vertices) != list:
            vertices = vertices.tolist()

        newMesh = Mesh.Mesh(vertices)
        return newMesh
    
    def isControlPoint(self, point, control_points):
        for control_point in control_points:
            if np.array_equal(point, control_point):
                return True
        return False
    
    #for attempt at control-point-based process
    def compute_weights(self, vertices, control_points):

        distances = np.linalg.norm(vertices[:, np.newaxis] - control_points, axis=2)
        inv_distances = np.reciprocal(distances, where=distances != 0)
        weights = inv_distances / np.sum(inv_distances, axis=1, keepdims=True)

        for i in range(len(vertices)): 
            #the points shouldn't move if the vertex's z coordinate is 0 (it's on the base), or if it's the control point
            #but at the same time we don't necessarily want every base-point to be a control point bc as shown that leads to weird behaviors
            #so we set the weight to 1.0 for some vertex's weights if we don't want the corresponding control point affecting it 
            if self.isControlPoint(vertices[i], control_points) or (vertices[i][2] < 0.1):
                print(f"Vertex shouldn't move: {vertices[i]}")
                weights[i] = 1.0
        
        #we also know we want the point closest to the top vertex that will exist, to stay near the top vertex
        #otherwise we won't have a point that stays there, which is how we ended up with the shape changing dimensions
        #so knowing the top vertex, we look for the closest point, and then set the weights for that corresponding point to 1.0
        top_vertex = np.array([5.0, 5.0, 10.0])
        closestToTop = np.argmin(np.linalg.norm(vertices - top_vertex, axis=1))
        weights[closestToTop] = 1.0

        return weights

    #use this one for now
    def pyramidFromCube(self, id='004'):

        print(f"Starting pyramid attempt")

        vertices = []

        for i in range(len(self.allmeshes[0].Facets)):
            facet_points = self.allmeshes[0].Facets[i].Points
            for point in facet_points:
                vertices.append(list(point))      

        vertices = np.array(vertices)

        # print(f"Original vertices: {vertices}")  
        # print(f"Number of vertices: {len(vertices)}")

        os.makedirs(f"pyramidtest{id}")
        meshUpdated = self.makeMesh(vertices)
        self.saveMeshSTL(meshUpdated, f"pyramidtest{id}/before_pyramid", self.meshres)

        faces = np.array([[facet.PointIndices[i] for i in range(3)] for facet in self.faces])
        print(f"Original faces: {faces}")    
        print(f"Number of faces: {faces.shape}")    

        control_points = np.array([
            [0, 10.0, 0],  #corner 3
            [10.0, 10.0, 0],   #corner 2
            [10.0, 0, 0],    #corner 1 
            [0, 0, 0],   #corner 0 
            [5.0, 5.0, 10.0], #top vertex
        ])

        weights = self.compute_weights(vertices, control_points)
        # print(f"Weights: {weights}")
        # print(f"Weights individual: {weights[0]}")

        count = 1

        deformed_vertices = vertices.copy()
        for i in range(len(vertices)):

            #shouldn't move the vertex if it's close enough to a control point
            #we previously set weights=1.0 for vertices where we don't want the corresponding control points to affect it
            if np.all(weights[i] != 1.0):
                deformed_vertices[i] = np.dot(weights[i], control_points)
                vertices[i] = np.dot(weights[i], control_points) 
                
            meshUpdated = self.makeMesh(deformed_vertices)
            mesh2 = self.makeMesh(vertices)

            if count % 500== 0: 
                self.saveMeshSTL(meshUpdated, f"pyramidtest{id}/pyramid_test_00{count/500}", self.meshres)
                self.saveMeshSTL(mesh2, f"pyramidtest{id}/pyramid_test_updatingvertex_00{count/500}", self.meshres)

            count += 1

        meshUpdated = self.makeMesh(deformed_vertices)
        print(f"Making new mesh: {meshUpdated}")

        print(f"New vertices: {deformed_vertices}")
        #print(f"New faces: {faces}")
        print(f"Number of new vertices: {len(deformed_vertices)}")
        #print(f"Number of new faces: {len(faces)}")

        self.saveMeshSTL(meshUpdated, f"pyramidtest{id}/pyramid_test_final", self.meshres)

        return deformed_vertices, faces

    def isPointOutsidePyramid(self, point, pyramidVertices):
        vertex0, vertex1, vertex2, vertex3, apex = pyramidVertices

        # Calculate face normals
        face0 = np.cross(vertex1 - vertex0, vertex2 - vertex0)
        face1 = np.cross(vertex2 - vertex1, vertex3 - vertex1)
        face2 = np.cross(vertex3 - vertex2, vertex0 - vertex2)
        face3 = np.cross(vertex0 - vertex3, vertex1 - vertex3)

        # Compute signed distances
        distances = [
            np.dot(point - vertex0, face0),
            np.dot(point - vertex1, face1),
            np.dot(point - vertex2, face2),
            np.dot(point - vertex3, face3)
        ]

        #if all distances are positive: point is inside pyramid
        #if any of the distances are negative: point is outside the pyramid
        if all(distance >= 0 for distance in distances):
            return True
        else:
            return False
        
    def pointDistanceToPyramid(self, point, pyramid_vertices):
        closest_dist = float("inf")
        closest_point = None

        for i in range(4):
            vertex0, vertex1, vertex2 = pyramid_vertices[i], pyramid_vertices[(i+1)%4], pyramid_vertices[4]

            #barycentric coordinates
            v0 = vertex2 - vertex0
            v1 = vertex1 - vertex0
            v2 = point - vertex0

            dot00 = np.dot(v0, v0)
            dot01 = np.dot(v0, v1)
            dot02 = np.dot(v0, v2)
            dot11 = np.dot(v1, v1)
            dot12 = np.dot(v1, v2)

            inv_denom = 1 / (dot00 * dot11 - dot01 * dot01)
            u = (dot11 * dot02 - dot01 * dot12) * inv_denom
            v = (dot00 * dot12 - dot01 * dot02) * inv_denom

            if u >= 0 and v >= 0 and u + v <= 1:
                # Point is inside the current face
                closest_point_on_face = vertex0 + u * v0 + v * v1
                dist = np.linalg.norm(point - closest_point_on_face)
                if dist < closest_dist:
                    closest_dist = dist
                    closest_point = closest_point_on_face

        return closest_dist, closest_point
    

    #alternative method: find boolean between mesh and target pyramid
    #then use that difference (the volume inside 1 but not the other, total) as objective fcn
    #need to write fcn for mesh (vertices) to solid in order to find boolean/intersection operation though for volume
    #which is annoying
    #so let's do v1 of that first? below
    
    def pyramidFromCubeV2(self, id='002'):

        print(f"Starting pyramid attempt")

        vertices = []

        for i in range(len(self.allmeshes[0].Facets)):
            facet_points = self.allmeshes[0].Facets[i].Points
            for point in facet_points:
                vertices.append(list(point))      

        vertices = np.array(vertices) #cube vertices now defined

        print(f"Original vertices: {vertices[0:3]}")  
        print(f"Number of vertices: {len(vertices)}")

        os.makedirs(f"pyramidtest{id}")
        meshUpdated = self.makeMesh(vertices)
        self.saveMeshSTL(meshUpdated, f"pyramidtest{id}/before_pyramid", self.meshres)

        faces = np.array([[facet.PointIndices[i] for i in range(3)] for facet in self.faces])
        # print(f"Original faces: {faces}")    
        print(f"Number of faces: {faces.shape}")   

        pyramidVertices = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 10.0, 0.0],
            [10.0, 10.0, 0.0],
            [10.0, 0.0, 0.0],
            [5.0, 5.0, 5.0],
        ])

        #v3 approach?
        #objective to minimize: (volume of mesh outside of target volume) + (volume of target outside of mesh), with intersect? 
        #determine pyramid volume based on defined vertices, & determine where this sits in space
        #how to calculate mesh volume outside of another volume? 
        #somehow move mesh points to reduce the volume: probably similar to grad desc? (how to determine though - move vertex to nearest point inside the mesh?)
        #eg: calculate the distance from a mesh vertex to the closest point inside the target area, then move it there (or just move to some distance within it)
        

        #for vertex in mesh: determine if the vertex is inside/outside target pyramid volume
        #either: minimize volume outside the target volume AND/OR: minimize number of vertex points outside the target volume (ie, want it to get to 0) - could be more straightforward?
        #keep moving the vertices until the number outside the pyramid target is 0

        #attempt1: iteratively moving vertices inside the pyramid if they aren't in there already
        count = 0
        for i in range(len(vertices)):
            vertex = vertices[i]
            if self.isPointOutsidePyramid(vertex, pyramidVertices):
                distance, targetPoint = self.pointDistanceToPyramid(vertex, pyramidVertices)
                print(f"Distance to pyramid surface: {distance}, target point: {targetPoint}")
                vertices[i] = targetPoint
            if count % 200 == 0: 
                meshUpdated = self.makeMesh(vertices)
                self.saveMeshSTL(meshUpdated, f"pyramidtest{id}/pyramid_test_00{count/200}", self.meshres)
            count += 1
        
        #attempt 2: calculating the intersection volume and using that

        # def objectiveFunction():
        #     return volumeDiff
        
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
    
    