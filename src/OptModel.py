import FreeCAD
import Part
import Mesh
import MeshPart
from FreeCAD import Base

class OptModel_Box:
    def __init__(self, threshold = 0.1, gprev = 10000, gcurr = 5, delstep = 0.05, del_e = 1):
        #all placeholder values for now I have no idea what they should be yet - maybe these shouldn't be here and instead in opt? 
        self.threshold_err = threshold
        self.g_prev = gprev
        self.g_curr = gcurr
        self.delstep = delstep
        self.del_e = del_e
        return
    
    def calculateDelE(self): 
        self.del_e = self.g_prev - self.g_curr
        return self.del_e

    
    def doTransform(self, cadModel, x=0, y=0, z=1): 

        print(f"delE: {self.del_e}, delStep: {self.delstep}, current g: {self.g_curr}")

        # del_theta = -1 * self.delstep * self.del_e #probably not the way

        del_theta = self.delstep * (self.g_curr - self.del_e)**2 * (-1 if (self.del_e > 0) else 1)

        del_theta_rot = FreeCAD.Rotation(0, z*del_theta, 0)

        axis = FreeCAD.Vector(x, y, z) 
        rot = FreeCAD.Placement(axis, del_theta_rot)
        # rot = FreeCAD.Rotation(axis, del_theta)
        print(f"Need to rotate by {del_theta}, transforming")

        cadModel.rotateByAmount(0, 0, del_theta, x, y, z)

        cadModel.CADdoc.recompute()

        return del_theta

    def updategValues(self, g_new):
        self.g_prev = self.g_curr
        self.g_curr = g_new
        return


class OptModel_3DRot:
    def __init__(self, threshold = 0.1, gprev = 10000, gcurr = 5, delstep = 0.1, del_e = 1):
        #all placeholder values for now I have no idea what they should be yet - maybe these shouldn't be here and instead in opt? 
        self.threshold_err = threshold
        self.g_prev = gprev
        self.g_curr = gcurr
        self.delstep = delstep
        self.del_e = del_e
        return

    def calculateError(self):
        return 

    # def calculateRotationMatrix(self):
    #     return 

    def calculateGradientRotation(self):
        return 

    def doTransform(self, cadModel, rotation):
        cadModel.rotateModel(rotation) #is rotation a rotation matrix? what exactly -> TBD -> rotation is a Matrix4D
        return


    # def gradientDescent3D(self, cadModel, rotationInit): #start off with whatever last one was, so maybe take in value?
    #     #set up loop... 
    #     return




    


class OptModel_Template:
    def __init__(self, threshold = 0.1, gprev = 10000, gcurr = 5, delstep = 0.1, del_e = 1):
        #all placeholder values for now I have no idea what they should be yet - maybe these shouldn't be here and instead in opt? 
        self.threshold_err = threshold
        self.g_prev = gprev
        self.g_curr = gcurr
        self.delstep = delstep
        self.del_e = del_e
        return
    
    def calculateDelE(self): 
        self.del_e = self.g_prev - self.g_curr
        return self.del_e

    
    def doTransform(self, cadModel, x=0, y=0, z=1): 
        return
    

    def updategValues(self, g_new):
        self.g_prev = self.g_curr
        self.g_curr = g_new
        return    

