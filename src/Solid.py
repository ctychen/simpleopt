import CADClass

class Box(CADClass.CAD): 
    
    # def __init__(self, stpfile):
    #     super(CADClass.CAD, self).__init__()
    #     self.STPfile = stpfile #want to be able to use as many HEAT cadclass fcns as possible, is this how to do it?
    #     return

    def __init__(self, stlfile, stpfile=""):
        super(CADClass.CAD, self).__init__()
        self.STLfile = stlfile
        self.STPfile = stpfile #want to be able to use as many HEAT cadclass fcns as possible, is this how to do it?
        return


    def createMesh(self): 
        #standard mesh for now, can change this
        meshes = self.load1Mesh(self.STLfile)
        #meshes = self.part2meshStandard(self.CAD)
        return meshes
    
