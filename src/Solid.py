import CADClass

class Solid(CADClass):
    def __init__(self, cadmodel):
        super(CADClass, self).__init__()
        self.cad = cadmodel


    def createMesh(self): 
        #standard mesh for now, can change this
        meshes = self.part2meshStandard(self.cad)
        return meshes
    
