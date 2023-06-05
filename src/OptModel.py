import FreeCAD
import Part
import Mesh
import MeshPart
from FreeCAD import Base
import numpy as np

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
    def __init__(self, objective, threshold = 0.1, gprev = 10000, gcurr = 5, delstep = 0.1, del_e = 1):
        #all placeholder values for now I have no idea what they should be yet - maybe these shouldn't be here and instead in opt? 
        self.objective = objective
        self.threshold_err = threshold
        self.g_prev = gprev
        self.g_curr = gcurr
        self.delstep = delstep
        self.del_e = del_e
        return

    def gradientDescent(self, cadModel, currentAng, calcQVal, angleRange):

        # velocity = momentum * velocity + learning_rate * gradient
        # next point = point - velocity

        # angleRange = 2 #could be user input

        # currentAng = cadModel.getCurrentRotationAngles() #this was from using the solid
        currentXAng = currentAng[0]
        currentYAng = currentAng[1]
        currentZAng = currentAng[2]

        xRot = np.linspace(currentXAng - angleRange, currentXAng + angleRange, angleRange*2)
        yRot = np.linspace(currentYAng - angleRange, currentYAng + angleRange, angleRange*2)
        zRot = np.linspace(currentZAng - angleRange, currentZAng + angleRange, angleRange*2)

        X, Y, Z = np.meshgrid(xRot, yRot, zRot)

        q_inGrid = np.zeros((angleRange*2, angleRange*2, angleRange*2))

        count = (angleRange*2)**3

        for i in range(len(xRot)):
                    for j in range(len(yRot)):
                        for k in range(len(zRot)):
                            xVal = xRot[i]
                            yVal = yRot[j]
                            zVal = zRot[k]
                            q_inGrid[i, j, k] = calcQVal(xVal, yVal, zVal)
                            print(f"Point done: {xVal}, {yVal}, {zVal}")
                            count -= 1
                            print(f"Points left: {count}")        

        q_1D = q_inGrid.flatten()

        min_q = np.amin(q_1D)
        idxMin = np.argmin(q_1D)

        xNew = X.flatten()[idxMin]
        yNew = Y.flatten()[idxMin]
        zNew = Z.flatten()[idxMin]

        return [[xNew, yNew, zNew], min_q]        


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

