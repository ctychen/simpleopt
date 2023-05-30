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
        return

    #todo: maybe can replace some stuff wiht objFromPartnum from CADClass directly

    def createMesh(self): 
        #standard mesh for now, can change this
        # meshes = self.load1Mesh(self.STLfile)
        self.parts = self.loadSTEP()
        self.meshes = self.part2meshStandard(self.CADparts)
        # meshSTL = self.writeMesh2file(meshes, "meshes")
        print(f"Meshed this:")
        print(self.meshes)
        print(type(self.meshes[0]))
        return self.meshes
    
