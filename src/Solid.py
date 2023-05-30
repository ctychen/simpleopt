import CADClass
import FreeCAD 
import Part

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

    def createMesh(self): 
        #standard mesh for now, can change this
        # meshes = self.load1Mesh(self.STLfile)
        # self.parts = self.loadSTEP()
        #self.meshes = self.part2meshStandard(self.CADparts)
        self.meshes = self.part2mesh(self.CADparts, 500)
        # meshSTL = self.writeMesh2file(meshes, "meshes")
        print("Meshed")
        # print(self.meshes)
        # print(type(self.meshes[0]))
        return self.meshes

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


        # axis = FreeCAD.Vector(x, y, z)#(1,0,0)
        # rot = FreeCAD.Rotation(thetax, thetay, thetaz)#(axis,ang)
        # self.CADdoc.Objects[0].Placement = FreeCAD.Placement(axis,rot)    
        # print("CAD Permutation Complete")
        return 
    
