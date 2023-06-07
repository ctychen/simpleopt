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

        #calculating the rotation matrix
        # rotationMatrix = np.eye(3)
        # rotationMatrix = np.dot(rotationMatrix, np.array([
        #     [1, 0, 0],
        #     [0, np.cos(xAngRad), -np.sin(xAngRad)],
        #     [0, np.sin(xAngRad), np.cos(xAngRad)]
        # ]))
        # rotationMatrix = np.dot(rotationMatrix, np.array([
        #     [np.cos(yAngRad), 0, np.sin(yAngRad)],
        #     [0, 1, 0],
        #     [-np.sin(yAngRad), 0, np.cos(yAngRad)]
        # ]))
        # rotationMatrix = np.dot(rotationMatrix, np.array([
        #     [np.cos(zAngRad), -np.sin(zAngRad), 0],
        #     [np.sin(zAngRad), np.cos(zAngRad), 0],
        #     [0, 0, 1]
        # ]))

        # print(f"Rotation matrix: {rotationMatrix}")

        #applying the rotation to the mesh coordinates
        # rotatedVertices = np.dot(vertices, rotationMatrix.T)   

        # print(f"Rotated vertices: {rotatedVertices}")   
        # 
        rotation = Rotation.from_euler('xyz', np.radians([xAng, yAng, zAng]), degrees=False)
        rotatedVertices = rotation.apply(vertices)

        return [rotation, rotatedVertices]
    
    def calculateRotationOnVector(self, vector, angles):
        #angles should be the target angle: rotate the vector to (angles), not by (angles)
        angles = np.radians(angles)
        r = Rotation.from_euler('xyz', angles)
        rotatedVector = r.apply(vector)        
        return rotatedVector
    
    # def calculateRotationOnVector_Fast(self, vector, angles):
    #     angles = np.radians(angles)

    #     theta_x, theta_y, theta_z = angles

    #     rotation_x = np.array([[1, 0, 0],
    #                         [0, np.cos(theta_x), -np.sin(theta_x)],
    #                         [0, np.sin(theta_x), np.cos(theta_x)]])

    #     rotation_y = np.array([[np.cos(theta_y), 0, np.sin(theta_y)],
    #                         [0, 1, 0],
    #                         [-np.sin(theta_y), 0, np.cos(theta_y)]])

    #     rotation_z = np.array([[np.cos(theta_z), -np.sin(theta_z), 0],
    #                         [np.sin(theta_z), np.cos(theta_z), 0],
    #                         [0, 0, 1]])

    #     rotatedVector = rotation_z @ (rotation_y @ (rotation_x @ vector))
    #     return rotatedVector
    
    def setVertices(self, newVertices):
        self.meshVertices = newVertices
        return 
    
    def findMeshCenter(self, mesh):
        center = mesh.BoundBox.Center
        return center
    
    
    def updateMesh(self, newVertices):
        # self.meshVertices = newVertices
        # # points = [FreeCAD.Vector(x, y, z) for x, y, z in newVertices]
        # # points = np.array(vertices)
        # newMesh = Mesh.Mesh(newVertices)
        # # newMesh.addPoints(points)
        # self.allmeshes[0] = newMesh

        if type(newVertices) != list:
            newVertices = newVertices.tolist()

        # self.verticesFromFacets = newVertices
        # points = [FreeCAD.Vector(x, y, z) for x, y, z in newVertices]
        # points = np.array(vertices)

        newMesh = Mesh.Mesh(newVertices)
        self.allmeshes[0] = newMesh
        return newMesh
    
    def makeMesh(self, vertices):
        # points = [FreeCAD.Vector(x, y, z) for x, y, z in vertices]
        # points = np.array(vertices)
        # print(type(vertices))
        # print(vertices.tolist())

        if type(vertices) != list:
            vertices = vertices.tolist()

        # print(len(triangles))
        newMesh = Mesh.Mesh(vertices)
        # newMesh.addPoints(points)
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
        # norms = self.calculateNorms(vertices)
        # norms = self.calculateNormals(vertices)
        # centers = self.calculateCenters(vertices)
        # areas = self.calculateAreas(vertices)    


        # for mesh in meshes:
        #     #mesh = obj.Mesh
        #     if (mesh == None) or (mesh=='None'):
        #         print("No Mesh for one of these objects.  Did you have a typo in input file?")
        #         print("Check HEAT output for Mesh Not Found errors")
        #         log.info("No Mesh for one of these objects.  Did you have a typo in input file?")
        #         log.info("Check HEAT output for Mesh Not Found errors")
        #     else:
        #         N_facets = mesh.CountFacets
        #         x = np.zeros((N_facets,3))
        #         y = np.zeros((N_facets,3))
        #         z = np.zeros((N_facets,3))

        #         for i,facet in enumerate(mesh.Facets):
        #             #mesh points
        #             for j in range(3):
        #                 x[i][j] = facet.Points[j][0]
        #                 y[i][j] = facet.Points[j][1]
        #                 z[i][j] = facet.Points[j][2]

        #         # scale and permute if necessary
        #         x,y,z = self.scale_and_permute(x,y,z)
        #         # get face normals and face centers
        #         norms.append(self.faceNormals(mesh))
        #         centers.append(self.faceCenters(x,y,z))
        #         areas.append(self.faceAreas(mesh))

        mesh = self.makeMesh(vertices)
        norms = self.faceNormals(mesh)
        centers = [] #self.faceCenters(mesh)
        areas = [] #self.faceAreas(mesh) 
        # print(f"Norms: {norms}")  #figure out how to process these things 
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
    

