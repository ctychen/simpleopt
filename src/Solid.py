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
        base = pyramidVertices[0:4]
        apex = pyramidVertices[4]
        
        for i in range(4):
            #vertices for pyramid edge
            v1, v2 = base[i], base[(i+1)%4]
            
            #find normal between edges
            normal = np.cross(v2 - v1, apex - v1)
            
            #vector from the point to the triangle face
            vec = v1 - point
            
            #if dot product is positive, the point is outside pyramid
            if np.dot(normal, vec) >= 0:
                return False
        
        return True
        
    def pointTriangleDistance(self, P, A, B, C):
        # Compute vectors        
        v0 = C - A
        v1 = B - A
        v2 = P - A

        # Compute dot products
        dot00 = np.dot(v0, v0)
        dot01 = np.dot(v0, v1)
        dot02 = np.dot(v0, v2)
        dot11 = np.dot(v1, v1)
        dot12 = np.dot(v1, v2)

        # Compute barycentric coordinates
        invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
        u = (dot11 * dot02 - dot01 * dot12) * invDenom
        v = (dot00 * dot12 - dot01 * dot02) * invDenom

        # Check if point is in triangle
        if (u >= 0) and (v >= 0) and (u + v <= 1):
            proj = A + u*v0 + v*v1
            return np.linalg.norm(P - proj), proj
        else:
            # Here, compute the distance to the edge or vertex
            edge_points = [
                (np.dot(P - A, B - A) / np.dot(B - A, B - A)) * (B - A) + A,
                (np.dot(P - B, C - B) / np.dot(C - B, C - B)) * (C - B) + B,
                (np.dot(P - C, A - C) / np.dot(A - C, A - C)) * (A - C) + C
            ]
            edge_dists = [np.linalg.norm(P - pt) for pt in edge_points]
            min_dist_idx = np.argmin(edge_dists)
            return edge_dists[min_dist_idx], edge_points[min_dist_idx]

    def pointDistanceToPyramid(self, P, pyramidVertices):
        base = pyramidVertices[:4]
        apex = pyramidVertices[4]
        
        results = [
            self.pointTriangleDistance(P, base[i], base[(i+1)%4], apex)
            for i in range(4)
        ]
        
        results.append(self.pointTriangleDistance(P, base[0], base[1], base[2]))
        results.append(self.pointTriangleDistance(P, base[2], base[3], base[0]))
        
        min_distance, closest_point = min(results, key=lambda x: x[0])
        
        return min_distance, closest_point
    
    def gradMoveVertex(self, vertex, pyramidVertices, step_size=0.01):
        distance, closest_point = self.pointDistanceToPyramid(vertex, pyramidVertices)
        #calculate direction vector from the vertex to the closest point on the pyramid
        #then normalize
        direction_vector = closest_point - vertex
        direction_vector /= np.linalg.norm(direction_vector)
        #move vertex towards closest point using direction vector
        vertex += step_size * direction_vector
        return vertex
    

    #alternative method: find boolean between mesh and target pyramid
    #then use that difference (the volume inside 1 but not the other, total) as objective fcn
    #need to write fcn for mesh (vertices) to solid in order to find boolean/intersection operation though for volume
    #which is annoying
    #so let's do v1 of that first? below
    
    def pyramidFromCubeV2(self, id='004'):

        print(f"Starting pyramid attempt")

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

        pyramidVertices = np.array([
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [10.0, 10.0, 0.0],
            [0.0, 10.0, 0.0],
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

        ####attempt1####: iteratively moving vertices inside the pyramid if they aren't in there already
        # count = 0
        # for i in range(len(vertices)):
        #     vertex = vertices[i]
        #     if self.isPointOutsidePyramid(vertex, pyramidVertices):
        #         distance, targetPoint = self.pointDistanceToPyramid(vertex, pyramidVertices)
        #         print(f"Distance to pyramid surface: {distance}, target point: {targetPoint}")
        #         vertices[i] = targetPoint
        #     if count % 200 == 0: 
        #         meshUpdated = self.makeMesh(vertices)
        #         self.saveMeshSTL(meshUpdated, f"pyramidtest{id}/pyramid_test_00{count/200}", self.meshres)
        #     count += 1
        
        #attempt 2: calculating the intersection volume and using that

        # def objectiveFunction(meshcube, pyramidMesh):
        #     make meshcube, pyramidMesh into solid
        #     do boolean with that, and calculate volume of boolean 
        #     return volumeDiff

        # def objectiveFunction(cubeVertices, pyramidVertices):
        #     cubeShape=Part.Shape()
        #     cubeShape.makeShapeFromMesh(self.makeMesh(cubeVertices),0.05) # the second arg is the tolerance for sewing
        #     cubeSolid = Part.makeSolid(cubeShape)
        #     pyramidShape=Part.Shape()
        #     pyramidShape.makeShapeFromMesh(self.makeMesh(pyramidVertices),0.05) # the second arg is the tolerance for sewing
        #     pyramidSolid = Part.makeSolid(pyramidShape)

        #     intersection = cubeSolid.Shape.common(pyramidSolid.Shape)
        #     volumeDiff = intersection.Volume

        #     return volumeDiff

        pVertices = [
            FreeCAD.Vector(0,0,0),
            FreeCAD.Vector(10,0,0),
            FreeCAD.Vector(10,10,0),
            FreeCAD.Vector(0,10,0),
            FreeCAD.Vector(5,5,10),  # the apex
        ]
        # Create the base
        pBase = Part.makePolygon(pVertices[:4] + [pVertices[0]])  # makePolygon requires a closed loop, hence append the first vertex at the end

        # Create the faces
        pFaces = [Part.Face(pBase)]
        for i in range(4):
            triangle = Part.makePolygon([pVertices[i], pVertices[(i+1)%4], pVertices[4], pVertices[i]])
            pFaces.append(Part.Face(triangle))

        # Create a shell from the faces
        shell = Part.makeShell(pFaces)

        # Create a solid from the shell
        pyramidSolid = Part.makeSolid(shell)
        
        psolid = FreeCAD.ActiveDocument.addObject("Part::Feature", "Pyramid")
        psolid.Shape = pyramidSolid

        # self.allmeshes = self.part2mesh(self.CADparts, meshres)

        # pyramidMesh = self.part2mesh(pyramidSolid, self.meshres)
        pyramidMesh = self.part2meshStandard(psolid)

        print(f'Pyramid Mesh vertices: {pyramidMesh[0].Points}')
        print(f'Pyramid Mesh faces: {pyramidMesh[0].Facets}')
        #above this: stuff works - created a pyramid mesh that can work

        #pyramidMeshVertices = pyramidMesh[0].Points
        #pyramidMeshFaces = pyramidMesh[0].Facets

        pyramidMeshVertices = []
        pyramidMeshFaces = []

        for i in range(len(pyramidMesh[0].Facets)):
            facet_points_list = pyramidMesh[0].Facets[i]
            facet_points = facet_points_list.Points
            for point in facet_points:
                pyramidMeshVertices.append(list(point))  
            for j in facet_points_list.Points:
                print(j)
                pyramidMeshFaces.append(list(j))

        pyramidMeshVertices = np.array(pyramidMeshVertices) 
        pyramidMeshFaces = np.array(pyramidMeshFaces)

        print(f"Pyramid mesh faces: {pyramidMeshFaces}")

        #below this: needs work 

        #in theory what this should do:
        #calculate the volume that's inside the cube/object but outside the target volume
        #calculate the volume that's inside the target voluem but outside the cube/object
        #sum the volume
        #this sum is what should be minimized
        #so it's time for mesh ops - boolean/intersection - and FC may not be the way to go

        # mesh = trimesh.Trimesh(vertices=[[0, 0, 0], [0, 0, 1], [0, 1, 0]],
        #                faces=[[0, 1, 2]])

        #both of these inputs should be trimesh objects
        #objective should be: (volume of mesh) - (volume of intersection of/boolean mesh & pyramid) -> minimize
        #if mesh fits perfectly inside pyramid, then Vmesh = Vpyr, so Vbool = Vmesh
        #so objective should be abs(Vmesh - Vbool(mesh,pyr)) since we want that diff to be close to 0
        #trimesh also has fcn for trimesh.boolean.difference
        #trimesh.boolean.intersection(meshes, engine=None, **kwargs)
        def objectiveFunction(cubeMesh, pyramidMesh):
            intersectionMesh = trimesh.boolean.intersection([cubeMesh, pyramidMesh])
            vIntersect = intersectionMesh.volume
            vMesh = cubeMesh.volume
            volumeDiff = np.abs(vMesh - vIntersect)
            return volumeDiff

        # def objectiveFunction(cubeVertices, pyramidSolid):
        #     # cubeShape=Part.Shape()
        #     # cubeMesh = self.makeMesh(cubeVertices)
        #     # cubeShape.makeShapeFromMesh(cubeMesh.Topology, 0.05) # the second arg is the tolerance for sewing
        #     # cubeSolid = Part.makeSolid(cubeShape)
            
        #     cubeMesh = self.makeMesh(cubeVertices)
        #     cubeShape = Part.makeSolid(Part.makeShell(cubeMesh.Facets))
        #     cubeSolid = FreeCAD.ActiveDocument.addObject("Part::Feature", "CubeSolid")
        #     cubeSolid.Shape = cubeShape

        #     intersection = cubeSolid.Shape.common(pyramidSolid.Shape)
        #     volumeMesh = cubeSolid.Shape.Volume
        #     volumeBool = intersection.Volume

        #     return volumeMesh - volumeBool
        
        cubeVertices = vertices.copy()
        cubeFaces = faces.copy()
        count = 0

        cubeTriMesh = trimesh.Trimesh(
            vertices=cubeVertices, faces=cubeFaces
        )

        pyramidTriMesh = trimesh.Trimesh(
            vertices=pyramidMeshVertices, faces=pyramidMeshFaces
        )

        pyramidTriMesh.export(f"pyramidtest{id}/basepyramid.stl")
        
        while objectiveFunction(cubeTriMesh, pyramidTriMesh) > 0:
            #move mesh vertices somehow
            #update the mesh somehow
            
            for i in range(len(cubeVertices)):
                #keep step size fixed for now but move vertex a bit, according to direction/distance from pyramid
                print(f"Original cube vertex {i}: {cubeVertices[i]}")

                cubeVertices[i] = self.gradMoveVertex(cubeVertices[i], pyramidVertices)
                cubeTriMesh.vertices[i] = cubeVertices[i] #Not sure if this is needed but have to check

                print(f"Modified cube vertex {i}: {cubeVertices[i]}")
                print(f"Modified cube vertex in Trimesh {i}: {cubeTriMesh.vertices[i]}")

                if count % 200 == 0:
                    #cubeMesh = self.makeMesh(cubeVertices)
                    #self.saveMeshSTL(cubeMesh, f"pyramidtest{id}/wip_{count}", self.meshres)
                    #mesh2.export('stuff.stl')
                    cubeTriMesh.export(f"pyramidtest{id}/wip_{count}.stl")

            count += 1


        #while objectiveFunction(meshcube, pyramidMesh) > 0:
        #   somehow move some mesh vertices 
        #   regenerate the meshes and set meshcube = new mesh after the modifications - distance change could be similar to control points?
        #   repeat process
        
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
    
    