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

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

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

    def gradientDescent(self, cadModel, currentAng, calcQVal, angleRange, plotEnable = True, savePlotsTo = '', ct = 0, numSamples = 24):

        currentXAng = currentAng[0]
        currentYAng = currentAng[1]
        currentZAng = currentAng[2]


        xRot = np.linspace(currentXAng - angleRange, currentXAng + angleRange, numSamples)
        yRot = np.linspace(currentYAng - angleRange, currentYAng + angleRange, numSamples)
        zRot = np.linspace(currentZAng - angleRange, currentZAng + angleRange, numSamples)

        q_inGrid = np.zeros((numSamples, numSamples, numSamples))

        count = (numSamples)**3

        for k in range(len(zRot)):
                    for j in range(len(yRot)):
                        for i in range(len(xRot)):
                            xVal = xRot[i]
                            yVal = yRot[j]
                            zVal = zRot[k]
                            q_inGrid[i, j, k] = calcQVal(xVal, yVal, zVal)
                            count -= 1            

        q_1D = q_inGrid.flatten()

        min_q = np.min(q_inGrid)
        min_indices = np.unravel_index(np.argmin(q_inGrid), q_inGrid.shape)
        xNew = xRot[min_indices[0]]
        yNew = yRot[min_indices[1]]
        zNew = zRot[min_indices[2]]
        

        if plotEnable:
            print(f"Min q value: {min_q}")
            fig = go.Figure(data = go.Surface(x = xRot, y = yRot, z = q_inGrid[:, :, min_indices[2]]))
            
            fig.update_layout(title_text=f"Min q value: {min_q}, at Z = {zRot[min_indices[2]]}")
            
            fig.show()            
            output_file = f"{savePlotsTo}/{ct}_step_qmin_{min_q}_at_z_{zRot[min_indices[2]]}.html"
            pio.write_html(fig, output_file)
            print(f"Plotted this iteration, saved file: {output_file}")

        return [[xNew, yNew, zNew], min_q, q_inGrid, q_1D, xRot, yRot, zRot]

                 


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

