import CADClass
import FreeCAD 
import Part
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

    def __init__(self, stlfile="", stpfile="", meshres=100):
        super(CADClass.CAD, self).__init__()
        self.STLfile = stlfile
        self.STPfile = stpfile
        self.parts = self.loadSTEP()
        # self.allmeshes = self.part2mesh(self.CADparts, meshres) #this returns a list of meshes - replaced self.meshes but name confusion
        # self.mesh = self.allmeshes[0] #this is a meshobject
        # self.meshVertices = np.array(self.mesh.Points)
        return
    

    def calculateRotationOnMesh(self, vertices, xAng, yAng, zAng):
        xAngRad = np.radians(xAng)
        yAngRad = np.radians(yAng)
        zAngRad = np.radians(zAng)

        #calculating the rotation matrix
        rotationMatrix = np.eye(3)
        rotationMatrix = np.dot(rotationMatrix, np.array([
            [1, 0, 0],
            [0, np.cos(xAngRad), -np.sin(xAngRad)],
            [0, np.sin(xAngRad), np.cos(xAngRad)]
        ]))
        rotationMatrix = np.dot(rotationMatrix, np.array([
            [np.cos(yAngRad), 0, np.sin(yAngRad)],
            [0, 1, 0],
            [-np.sin(yAngRad), 0, np.cos(yAngRad)]
        ]))
        rotationMatrix = np.dot(rotationMatrix, np.array([
            [np.cos(zAngRad), -np.sin(zAngRad), 0],
            [np.sin(zAngRad), np.cos(zAngRad), 0],
            [0, 0, 1]
        ]))

        #applying the rotation to the mesh coordinates
        rotatedVertices = np.dot(vertices, rotationMatrix.T)        

        return [rotationMatrix, rotatedVertices]
    #should we return the calculated rotation matrix, or the resulting list of points once the thing is applied? or should this do the rotation - maybe not

    #calculating norm, center, area if we start with mesh vertices
    def calculateNorms(self, vertices):
        return np.linalg.norm(vertices, axis=1)
    
    def calculateCenters(self, vertices):
        return np.mean(vertices, axis=1)
    
    def calculateAreas(self, vertices):
        v0 = vertices[:, 0, :]
        v1 = vertices[:, 1, :]
        v2 = vertices[:, 2, :]
        return 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)

    def setVertices(self, newVertices):
        self.meshVertices = newVertices
        return 
    
    def updateMesh(self, newVertices):
        points = [FreeCAD.ActiveDocument.Vector(x, y, z) for x, y, z in newVertices]
        newMesh = FreeCAD.Mesh.Mesh()
        newMesh.addPoints(points)
        return
    
    def processModel(self): 
        #no need to re-mesh if we're doing vertices 
        self.norms = self.calculateNorms(self.meshVertices)
        self.centers = self.calculateCenters(self.meshVertices)
        self.areas = self.calculateAreas(self.meshVertices)
        print("Found norms, centers, areas for vectorized setup")
        return

    # def updateMeshVertices(self, newVertices):
    #     # self.
    #     return

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
    
