import time
import sys
import numpy as np
import os

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots


import matplotlib.pyplot as plt

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
        self.opt = OptModel.OptModel_3DRot(g_obj)

        self.q_max_threshold = 9 #set this to whatever

        #self.box.loadSTEP()
        # mesh = self.box.load1Mesh(stlPath)

        self.Nx = 3
        self.Ny = 3
        self.Nz = 3

        self.xRot = np.linspace(-45.0, 45.0, self.Nx)
        self.yRot = np.linspace(-45.0, 45.0, self.Ny)
        self.zRot = np.linspace(-45.0, 45.0, self.Nz)

        self.xAng, self.yAng, self.zAng = np.meshgrid(self.xRot, self.yRot, self.zRot)

        self.del_theta = 0
        return

    def calcPeakQWithRotation(self, xR, yR, zR):
        self.box.rotateTo(xR, yR, zR) #use this rotate fcn
        self.fwd.processCADModel()
        q_mesh_all = self.fwd.calcQMesh()
        qPeak = max(q_mesh_all)
        print(f"Calculated Peak Q: {qPeak}")
        return qPeak


    def runModel(self, stlEn = True, plotEn = True):

        all_q_found = [] #i-th index will correspond to qval calculated on i-th iteration, before the transform
        all_rotations_found = [] #i-th index will correspond to rotation found on i-th iteration, so the i-th transform
        noVals = True

        while (noVals or (np.amin(all_q_found) > self.q_max_threshold)): 

            if noVals: noVals = False

            # self.fwd.processCADModel()

            # q_mesh_all = self.fwd.calcQMesh()

            # g_now = self.fwd.calcObjective(q_mesh_all)
            # print(f"g value found: {g_now}")

            # self.opt.updategValues(g_now)

            resultFromStep = self.opt.gradientDescent(self.box, self.calcPeakQWithRotation)
            #apply this rotation and repeat
            min_q_result = resultFromStep[1]
            rotation_result = resultFromStep[0] #format: [x, y, z]

            self.box.rotateTo(rotation_result[0], rotation_result[1], rotation_result[2])
            self.fwd.processCADModel()

            all_q_found.append(min_q_result)
            all_rotations_found.append(rotation_result)

        #once that condition reached, should mean that we're done optimizing, and so we can export a final file
        self.box.CADdoc.recompute()
        self.box.saveSTEP("final_box_3drot.step", self.box.CADobjs)
        
        return


    def plotRotations(self):

        print(f"Total number of points: {len(self.xAng) * len(self.yAng) * len(self.zAng)}")
        total = len(self.xAng) * len(self.yAng) * len(self.zAng)
        count = total
        qPeak_all = np.zeros((self.Nx, self.Ny, self.Nz))

        for i in range(len(self.xRot)):
            for j in range(len(self.yRot)):
                for k in range(len(self.zRot)):
                    xVal = self.xRot[i]
                    yVal = self.yRot[j]
                    zVal = self.zRot[k]
                    newQPeak = self.calcPeakQWithRotation(xVal, yVal, zVal)
                    qPeak_all[i, j, k] = self.calcPeakQWithRotation(xVal, yVal, zVal)
                    print(f"Point done: {xVal}, {yVal}, {zVal}")
                    count -= 1
                    print(f"Points left: {count}")


        qPeak_1D = qPeak_all.flatten()
        globalMinQ = np.amin(qPeak_1D) 
        idxMin = np.argmin(qPeak_1D)
        print(f"Global minimum of max(q): {globalMinQ} at index: {idxMin}")

        xFlat=self.xAng.flatten(),
        yFlat=self.yAng.flatten(),
        zFlat=self.zAng.flatten(),

        print(qPeak_all)

        xMin= xFlat[0][idxMin]
        yMin= yFlat[0][idxMin]
        zMin= zFlat[0][idxMin]

        print(f"Angles for minimum: {xMin}, {yMin}, {zMin}")

        # Create a 3D surface plot with color mapped to function values
        fig_4d = go.Figure(data=go.Volume(
            x=self.xAng.flatten(),
            y=self.yAng.flatten(),
            z=self.zAng.flatten(),
            value=qPeak_all.flatten(),
            isomin=np.min(qPeak_all),
            isomax=np.max(qPeak_all),
            opacity=0.3,
            surface_count=17,
            colorscale='Thermal'
        ))

        # # Set plot layout and axis labels
        fig_4d.update_layout(
            title='Peak Heat Flux Across Rotations',
            scene=dict(
                xaxis_title='X',
                yaxis_title='Y',
                zaxis_title='Z'
            )
        )

        # zIndex = np.abs(self.zRot - zMin).argmin()
        # q_zFixed = qPeak_all[:, :, zIndex]
        # fig_z = go.Figure(data=go.Volume(
        #     x=self.xAng.flatten(),
        #     y=self.yAng.flatten(),
        #     z=(q_zFixed.flatten(),),
        #     value=q_zFixed.flatten(),
        #     isomin=np.min(q_zFixed),
        #     isomax=np.max(q_zFixed),
        #     opacity=0.3,
        #     surface_count=17,
        #     colorscale='Thermal'
        # ))

        # # Set plot layout and axis labels
        # fig_z.update_layout(
        #     title='Peak Heat Flux Across Rotations',
        #     scene=dict(
        #         xaxis_title='X',
        #         yaxis_title='Y',
        #         zaxis_title='Z'
        #     )
        # )

        # Show the plot
        fig_4d.show()
        # fig_z.show()

        output_file_4d = 'plot_attempt_4d.html'
        pio.write_html(fig_4d, output_file_4d)

        # pio.write_html(fig_z, 'plot_z.html')

        print(f"Plotted Rotations Space")
        return globalMinQ


        

if __name__ == '__main__':

    t0 = time.time()

    # setup = RunSetup_1DBox()
    setup = RunSetup_3DBox()
    setup.runModel()
    # setup.plotRotations()

    print(f"Time elapsed: {time.time() - t0}")





        


    
