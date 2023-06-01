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

    def __init__(self, stlfile, stpfile=""):
        super(CADClass.CAD, self).__init__()
        self.STLfile = stlfile
        self.STPfile = stpfile
        self.parts = self.loadSTEP()
        return

    #todo: maybe can replace some stuff wiht objFromPartnum from CADClass directly

    def createMesh(self, res=1000): 
        self.meshes = self.part2mesh(self.CADparts, res)
        print("Meshed")
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
    
