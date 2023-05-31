import time
import sys
import numpy as np
import os

import plotly.graph_objects as go
import plotly.io as pio

try:
    runMode = os.environ["runMode"]
except:
    runMode = 'local'
    os.environ["runMode"] = runMode

if (runMode == 'docker'):
    FreeCADPath = '/usr/lib/freecad-daily/lib'
    # HEATPath = '/root/source/HEAT'
else:
    FreeCADPath = '/usr/lib/freecad-python3/lib'
    # HEATPath = '/Users/cchen/Desktop/HEAT'

sys.path.append(FreeCADPath)
print(sys.path)

import CADClass

import Solid
import ForwardModel
import OptModel

class RunSetup_1DBox:
    def __init__(self):
        g_obj = lambda qvals: max(qvals) #+ qvals.count(max(qvals)) #maybe changing obj function helps??

        stpPath = "test_box.step" 
        stlPath = " " #"box.stl"
        qDirIn = [0.0, -1.0, 0] #[m]
        qMagIn = 10.0 #[W/m^2]

        self.box = Solid.Box(stlPath, stpPath) 
        self.fwd = ForwardModel.ForwardModel_Box(g_obj, self.box, qMagIn, qDirIn) 
        self.opt = OptModel.OptModel_Box()

        #self.box.loadSTEP()
        # mesh = self.box.load1Mesh(stlPath)

        self.del_theta = 0
        return

    def runModel(self, stlEn = True):
        count = 0
        while ((abs(self.opt.del_e) > self.opt.threshold_err)):
            self.fwd.processCADModel()
            q_mesh_all = self.fwd.calcQMesh()
            g_now = self.fwd.calcObjective(q_mesh_all)
            print(f"g value found: {g_now}")
            self.opt.updategValues(g_now)
            self.opt.calculateDelE()

            if (abs(self.opt.del_e) < self.opt.threshold_err) and (self.opt.g_prev > self.opt.g_curr):
                print(f"[err]: {self.opt.del_e} [g_x now]: {self.opt.g_curr} [g_x prev]: {self.opt.g_prev} [theta]: {self.del_theta}")
                print(f"found opt, last err: {self.opt.del_e}, rotated: {self.del_theta}")
                self.box.CADdoc.recompute()
                self.box.saveSTEP("final_box.step", self.box.CADobjs)
                break
            else: 
                print(f"transform needed, error: {self.opt.del_e}")
                self.del_theta = self.opt.doTransform(self.box) #this prob doesn't match yet, gonna fix
                print(f"transformed: [err]: {self.opt.del_e} [g_x now]: {self.opt.g_curr} [g_x prev]: {self.opt.g_prev} [theta]: {self.del_theta}")
                
                if stlEn:
                    self.box.saveMeshSTL(self.box.meshes, f"boxmesh_{count}", 500)

            count += 1
        print(f"Iterations: {count}")
        return


class RunSetup_3DBox:
    def __init__(self):
        g_obj = lambda qvals: max(qvals) #+ qvals.count(max(qvals)) #maybe changing obj function helps??

        stpPath = "test_box.step" 
        stlPath = " " #"box.stl"
        qDirIn = [0.0, -1.0, 0] #[m]
        qMagIn = 10.0 #[W/m^2]

        self.box = Solid.Box(stlPath, stpPath) 
        self.fwd = ForwardModel.ForwardModel_Box(g_obj, self.box, qMagIn, qDirIn) 
        self.opt = OptModel.OptModel_3DRot()

        self.q_max = 9 #set this to whatever

        #self.box.loadSTEP()
        # mesh = self.box.load1Mesh(stlPath)

        self.Nx = 10
        self.Ny = 10
        self.Nz = 10

        self.xRot = np.linspace(-45.0, 45.0, self.Nx)
        self.yRot = np.linspace(-45.0, 45.0, self.Ny)
        self.zRot = np.linspace(-45.0, 45.0, self.Nz)

        self.xAng, self.yAng, self.zAng = np.meshgrid(self.xRot, self.yRot, self.zRot)

        self.del_theta = 0
        return

    def gradientDescent3D(self, iterations):
        
        # for i in range(iterations):
        #     flux_mesh_all = self.fwd.calcQMesh()
        #     gradRot = self.opt.calculateGradientRotation() #returns list with both


        return 

    def calcPeakQWithRotation(self, xR, yR, zR):
        # self.box.rotateByAmount()
        # initialMatrix = Matrix4D() #starts with an identity matrix since not sure what else to do
        #syntax is Rotation(rotateFrom, rotateTo)

        # for obj in FreeCAD.ActiveDocument.Objects:
        #     if type(obj) == Part.Feature:

                # print(f"Before Modifying Placement: {obj.Placement}")
                # obj.Placement = rot.multiply(obj.Placement)
                # obj.recompute()
                # print(f"After Modifying Placement: {obj.Placement}")
        # rotMatrix = FreeCAD.Rotation(FreeCAD.Vector(x, y, z), FreeCAD.Vector(xR, yR, zR))
        # currentPos = self.box.rotateModel()

        #make new mesh after applying the rotation, which is from current rotation to newly defined one. All coordinates are angles - Euler angles?
        #calculate q_all for this mesh
        #find max of q_all

        #how to make sure we start in the same place each time? Should all be an absolute rotation from (0,0,0) and the same initial fixed rotation.... 
        #maybe in this case we do have to re-load-in the model each time? 

        self.box.rotateTo(xR, yR, zR)
        self.fwd.processCADModel()
        q_mesh_all = self.fwd.calcQMesh()
        qPeak = max(q_mesh_all)
        print(f"Calculated Peak Q: {qPeak}")
        return qPeak

    def plotRotations(self):

        print(f"Total number of points: {len(self.xAng) * len(self.yAng) * len(self.zAng)}")
        count = len(self.xAng) * len(self.yAng) * len(self.zAng)
        qPeak_all = np.zeros((self.Nx, self.Ny, self.Nz))
        # print(self.xRot)
        # print(self.yRot)
        # print(self.zRot)

        # for i in len(self.xRot):
        #     for j in len(self.yRot):
        #         for k in len(self.zRot):
        #             qPeak_all[i][j][k] = self.calcPeakQWithRotation() 
        #             print(qPeak_all[i][j][k])

        for i in range(len(self.xAng)):
            for j in range(len(self.yAng)):
                for k in range(len(self.zAng)):
                    pointX = self.xAng[i, j, k]
                    pointY = self.yAng[i, j, k]
                    pointZ = self.zAng[i, j, k]
                    pointX2 = self.xRot[i]
                    pointY2 = self.yRot[j]
                    pointZ2 = self.zRot[k]
                    print(f"{pointX}, {pointY}, {pointZ}")
                    print({f"{pointX2}, {pointY2}, {pointZ2}"})
                    qPeak_all[i, j, k] = self.calcPeakQWithRotation(pointX, pointY, pointZ)
                    # print(qPeak_all[i, j, k])
                    count -= 1
                    print(f"Points left: {count}")

        fig = go.Figure()
        #fig.add_trace(go.Scatter3d(self.xRot, self.yRot, self.zRot, qPeak_all))

        # fig = go.Figure(data=go.Surface(x=self.xRot, y=self.yRot, z=self.zRot, surfacecolor=qPeak_all))

        # surface = go.Surface(x=self.xRot.flatten(), y=self.yRot.flatten(), z=self.zRot.flatten(), surfacecolor=qPeak_all.flatten())
        # contour = go.Contour(x=self.xRot.flatten(), y=self.yRot.flatten(), z=self.zRot.flatten(), value=qPeak_all.flatten())

        # fig = go.Figure(data=[surface])

        fig.add_trace(go.Surface(x=self.xAng.flatten(), y=self.yAng.flatten(), z=self.zAng.flatten(), surfacecolor=qPeak_all))

        # Set plot layout and axis labels
        fig.update_layout(
            title='Function Values with Contours in 3D Grid',
            scene=dict(
                xaxis_title='X Angle',
                yaxis_title='Y Angle',
                zaxis_title='Z Angle',
                xaxis = dict(nticks=4, range=[-45.0, 45.0],),
                yaxis = dict(nticks=4, range=[-45.0, 45.0],),
                zaxis = dict(nticks=4, range=[-45.0, 45.0],),
            )
        )

        # fig.update_layout(xaxis_range[-45.0, 45.0])
        # fig.update_layout(yaxis_range[-45.0, 45.0])
        # fig.update_layout(zaxis_range[-45.0, 45.0])

        # Show the plot
        fig.show()

        output_file = 'plot.html'
        pio.write_html(fig, output_file)

        print(f"Plotted Rotations Space")
        return

    


    def runModel(self, stlEn = True, plotEn = True):

        all_q_max = [] #i-th index will correspond to qval calculated on i-th iteration, before the transform
        all_rotations_found = [] #i-th index will correspond to rotation found on i-th iteration, so the i-th transform

        while (max(all_q_max) > self.q_max_threshold): 

            self.fwd.processCADModel()

            q_mesh_all = self.fwd.calcQMesh()

            all_q_max.append(max(q_mesh_all)) 

            g_now = self.fwd.calcObjective(q_mesh_all)
            print(f"g value found: {g_now}")

            self.opt.updategValues(g_now)

            rotation_found = self.gradientDescent3D() #to be written, but calculate grad desc - this will have an internal loop that runs x times

            all_rotations_found.append(rotation_found)

        #once that condition reached, should mean that we're done optimizing, and so we can export a final file
        self.box.CADdoc.recompute()
        self.box.saveSTEP("final_box_3drot.step", self.box.CADobjs)
        
        return

            # if (abs(self.opt.del_e) < self.opt.threshold_err) and (self.opt.g_prev > self.opt.g_curr):
            #     print(f"[err]: {self.opt.del_e} [g_x now]: {self.opt.g_curr} [g_x prev]: {self.opt.g_prev} [theta]: {self.del_theta}")
            #     print(f"found opt, last err: {self.opt.del_e}, rotated: {self.del_theta}")
            #     self.box.CADdoc.recompute()
            #     self.box.saveSTEP("final_box.step", self.box.CADobjs)
            #     break
            # else: 
            #     print(f"transform needed, error: {self.opt.del_e}")
            #     self.del_theta = self.opt.doTransform(self.box) #this prob doesn't match yet, gonna fix
            #     print(f"transformed: [err]: {self.opt.del_e} [g_x now]: {self.opt.g_curr} [g_x prev]: {self.opt.g_prev} [theta]: {self.del_theta}")
                      
        

if __name__ == '__main__':

    t0 = time.time()

    # setup = RunSetup_1DBox()
    setup = RunSetup_3DBox()
    # setup.runModel()
    setup.plotRotations()

    print(f"Time elapsed: {time.time() - t0}")





        


    
