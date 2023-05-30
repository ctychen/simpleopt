import FreeCAD
import Part
import Mesh
import MeshPart
from FreeCAD import Base

class OptModel_Box:
    def __init__(self, threshold = 0.5, gprev = 10000, gcurr = 1):
        #all placeholder values for now I have no idea what they should be yet - maybe these shouldn't be here and instead in opt? 
        self.threshold_err = 0.5
        self.g_prev = gprev
        self.g_curr = gcurr
        self.delstep = 5
        self.del_e = 1
        return
    
    def calculateDelE(self): 
        self.del_e = self.g_prev - self.g_curr
        return self.del_e
    
    def doTransform(self, cadModel, x=0, y=0, z=1): 
        # theta_new = theta - delstep*del_e -> angle input for rot is how much to rotate by so just delstep*del_e
        #calculate angle to rotate by
        #define vector for axis to rotate about: could put some constraints (or maybe we take than as an input? what's better for modularity?)
        #create freecad rotation object - what this needs is a vector [x,y,z] with x/y/z=1 for rotate about this axis, 0 for don't, and an angle
        #either that or use euler angles - eventually maybe, or just have 3 separate rotation objects, 1 for each axis, for more individual control?
        #keep box in same place as constraint: could use placement?
        #do the rotation (permutation on the object)
        #after this, cad file will have been modified

        del_theta = -1 * self.delstep * self.del_e
        axis = FreeCAD.Vector(x, y, z) 
        rot = FreeCAD.Rotation(axis, del_theta)
        print(f"Need to rotate by {del_theta}, transforming")

        # center = FreeCAD.Vector(0, 0, 0) #assuming cube centered on origin, doesn't really matter for simple rotation but maybe later
        #looks like for Placement you can define a local origin which would be what center would be - 
        #to do that you add pos = Placement.Base and then add pos to below 
        transform = FreeCAD.Placement(axis, rot)

        transform_base = Base.Rotation(axis, del_theta)

        # cadModel.CADdoc.recompute() #needed here?
        #for obj in cadModel.Objects: #to be fixed, take in mesh
        # print(cadModel.meshes[0])

        #attempt 1...

        # for obj in cadModel.meshes:
        #     print(obj)
        #     obj.Placement = transform.multiply(obj.Placement)
        #     if type(obj) == Part.Feature:
        #         obj.Placement = transform.multiply(obj.Placement)
        #         #cadModel.CADdoc.recompute() #not sure if we need this?
        #         #print(f"Recomputed with transformation on mesh")
        #         print(f"Applied transform?") #seems like code isn't getting here, wack
        #     print(f"Applied transform for object")
        # cadModel.CADdoc.recompute()

        #attempt 2... (maybe this is not how mesh elements work?)

        # print(dir(cadModel.meshes[0]))

        # face_list = cadModel.meshes[0].Facets
        # for face in face_list:
        #     # face.Placement = transform.multiply(face.Placement)
        #     face.rotate(axis, del_theta) #or, could do this?
        #     print(f"Applied transform")
        #     cadModel.CADdoc.recompute()

        #attempt 3...

        mesh = cadModel.meshes[0] 

        mesh.Placement.Rotation = mesh.Placement.Rotation.multiply(transform_base)

        print("Transformed?????")

        cadModel.CADdoc.recompute()

        return del_theta
    

    def updategValues(self, g_new):
        self.g_prev = self.g_curr
        self.g_curr = g_new
        return

