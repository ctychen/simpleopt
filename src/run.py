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

        self.box = Solid.Box_Vector(stlPath, stpPath) 
        self.fwd = ForwardModel.ForwardModel_Box(g_obj, self.box, qMagIn, qDirIn) 
        self.opt = OptModel.OptModel_3DRot(g_obj)

        # self.q_max_threshold = 6.0 #set this to whatever, 6 might be ambitious??

        #self.box.loadSTEP()
        # mesh = self.box.load1Mesh(stlPath)

        self.Nx = 100
        self.Ny = 100
        self.Nz = 100

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
        # print(f"Calculated Peak Q: {qPeak}")
        return qPeak
    
    def calcPeakQWithRotation_Vector(self, xR, yR, zR):
        #vector rotation math
        #process vectors/mesh to get the norms, centers, areas
        #calc mesh from there with lists of each
        #returned: [rotationMatrix, rotatedVertices]
        rotationResults = self.box.calculateRotationOnMesh(self.box.meshVertices, xR, yR, zR)
        #rotationMatrix = rotationResults[0]
        rotatedVertices = rotationResults[1]
        #self.fwd.process() -> technically don't need this here since not updating the mesh until we find a minimum
        q_mesh_all = self.fwd.calcQMesh_Vector(rotatedVertices)
        qPeak = max(q_mesh_all)
        return qPeak


    def runModel(self, momentum, epsilon_q = 0.01, epsilon_vel = 0.01, angleRange = 2, startingAngles = [0, 0, 0], stlEn = True, plotEn = True):

        all_q_found = [] #i-th index will correspond to qval calculated on i-th iteration, before the transform
        all_rotations_found = [] #i-th index will correspond to rotation found on i-th iteration, so the i-th transform
        noVals = True
        count = 0
        prev_q = 1000
        curr_q = 0
        prev_vel = 0
        curr_vel = 0
        angles = startingAngles
        # below_margin_count = 0

        #todo add momentum and maybe add variable learning rate
        # poss new condition: while del_velocity > threshold and del_q > threshold: run gradient descent. once both conditions met, done
        #while (noVals or (np.amin(all_q_found) > threshold)): 
        while (noVals or ((abs(prev_q - curr_q) > epsilon_q) and (abs(prev_vel - curr_vel) > epsilon_vel))):

            #picked something arbitrary for angleRange but idea would be enlarge area we're looking at each step, not sure if this is the way to do it? 
            #if difference between successive is large, take bigger steps/look at bigger area? should i even use delstep here? 

            resultFromStep = self.opt.gradientDescent(self.box, angles, self.calcPeakQWithRotation_Vector, curr_vel * self.opt.delstep) #eventually, this could be something completely different - look into pymoo

            #apply this rotation and repeat
            min_q_result = resultFromStep[1]
            rotation_result = resultFromStep[0] #format: [x, y, z]

            angles = rotation_result #since we can't get Rotation object easily from list of vertices

            prev_q = curr_q
            curr_q = min_q_result

            prev_vel = curr_vel 
            curr_vel = momentum * curr_vel + self.opt.delstep * (prev_q - curr_q)

            #save state of box before any rotations applied, so on 1st time

            if stlEn and noVals:
                #self.box.saveMeshSTL(self.box.meshes, f"outputs/boxmesh_initial", 500)
                self.box.saveMeshSTL(self.box.allmeshes, f"outputs/boxmesh_initial", 500)

            #result from graddesc will be list of rotated vertices
            #generate mesh from resulting vertices
            #feed this into next 

            #self.box.rotateTo(rotation_result[0], rotation_result[1], rotation_result[2])
            #self.fwd.processCADModel()

            #apply rotation and save updated model 
            # [rotationMatrix, rotatedVertices] returned from  def calculateRotationOnMesh(self, vertices, xAng, yAng, zAng):

            rotationResults = self.box.calculateRotationOnMesh(self.box.meshVertices, angles[0], angles[1], angles[2])
            verticesRotated = rotationResults[1]
            self.box.updateMesh(verticesRotated)
            
            if stlEn:
                #self.box.saveMeshSTL(self.box.meshes, f"outputs/boxmesh_{count}", 500)
                self.box.saveMeshSTL(self.box.allmeshes, f"outputs/{momentum}_run/boxmesh_{count}", 500)

            all_q_found.append(min_q_result)
            all_rotations_found.append(rotation_result)

            if noVals: noVals = False

            count += 1

        #once that condition reached, should mean that we're done optimizing, and so we can export a final file
        self.box.CADdoc.recompute()
        
        if stlEn:
                #self.box.saveMeshSTL(self.box.meshes, f"outputs/boxmesh_final", 500)
                self.box.saveMeshSTL(self.box.allmeshes, f"outputs/boxmesh_final", 500)

        self.box.saveSTEP("outputs/final_box_3drot.step", self.box.CADobjs)

        print(f"Reached below threshold in {count} iterations, minimum q is {all_q_found[len(all_q_found) - 1]}")
        print(f"Initial rotation: {all_rotations_found[0]}, found best rotation: {all_rotations_found[len(all_rotations_found) - 1]}")
        
        return [all_q_found, all_rotations_found]


    def plotRotations(self):

        # print(f"Total number of points: {len(self.xAng) * len(self.yAng) * len(self.zAng)}")
        total = len(self.xAng) * len(self.yAng) * len(self.zAng)
        count = total
        qPeak_all = np.zeros((self.Nx, self.Ny, self.Nz))

        for i in range(len(self.xRot)):
            for j in range(len(self.yRot)):
                for k in range(len(self.zRot)):
                    xVal = self.xRot[i]
                    yVal = self.yRot[j]
                    zVal = self.zRot[k]
                    # newQPeak = self.calcPeakQWithRotation(xVal, yVal, zVal)
                    qPeak_all[i, j, k] = self.calcPeakQWithRotation(xVal, yVal, zVal)
                    # print(f"Point done: {xVal}, {yVal}, {zVal}")
                    count -= 1
                    # print(f"Points left: {count}")


        qPeak_1D = qPeak_all.flatten()
        globalMinQ = np.amin(qPeak_1D) 
        idxMin = np.argmin(qPeak_1D)
        print(f"Global minimum of max(q): {globalMinQ} at index: {idxMin}")

        xFlat=self.xAng.flatten(),
        yFlat=self.yAng.flatten(),
        zFlat=self.zAng.flatten(),

        # print(qPeak_all)

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

        output_file_4d = 'plot_attempt_4d_2.html'
        pio.write_html(fig_4d, output_file_4d)

        # pio.write_html(fig_z, 'plot_z.html')

        print(f"Plotted Rotations Space")
        return globalMinQ


        

if __name__ == '__main__':

    t0 = time.time()

    # setup = RunSetup_1DBox()
    setup = RunSetup_3DBox()
    # all_q_found = setup.runModel(threshold=5.88)
    setup.plotRotations()

    # print(f"Initial heat flux: {all_q_found[0][0]}, best heat flux: {min(all_q_found[0])}")

    print(f"Time elapsed: {time.time() - t0}")





        


    
