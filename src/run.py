import time
import sys
import numpy as np
import os

import pandas as pd

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import plotly.express as px

import matplotlib.pyplot as plt


try:
    runMode = os.environ["runMode"]
except:
    runMode = 'local'
    os.environ["runMode"] = runMode

if (runMode == 'docker'):
    FreeCADPath = '/usr/lib/freecad-daily/lib'
    BlenderPath = '/usr/lib/blender'
    OpenSCADPath = '/usr/bin/openscad'

    # HEATPath = '/root/source/HEAT'
else:
    #FreeCADPath = '/Applications/FreeCAD.app' #'/usr/lib/freecad-python3/lib'
    FreeCADPath = '/Applications/FreeCAD.app/Contents/Resources/lib'
    BlenderPath = '/Applications/Blender.app/Contents/MacOS/blender'
    # HEATPath = '/Users/cchen/Desktop/HEAT'

sys.path.append(FreeCADPath)
sys.path.append(BlenderPath)
sys.path.append(OpenSCADPath)
print(sys.path)

import trimesh

import CADClass

import Solid
import ForwardModel
import OptModel

class RunSetup_MeshHF:
    def __init__(self):
        g_obj = lambda qvals: max(qvals) #+ qvals.count(max(qvals)) #maybe changing obj function helps??

        stpPath = "unit_test_cube.step" #"unit_test_cone.step" 
        # stpPath = "test_sphere.step"
        # stpPath = "test_pfc_block.step" #"unit_test_pfc.step" #for multiple directions 

        stlPath = " " #"box.stl"

        # qDirIn = [[0.707, -0.707, 0.0], [-0.707, -0.707, 0.0]] #[m]
        qDirIn = [0.0, -1.0, 0.0] #[m]
        # qDirIn = [0.707, -0.707, 0.0] #[m]
        # qDirIn = [-0.707, -0.707, 0.0] #[m]
        # qDirIn = [0.0, -0.707, 0.707] #[m]
        # qDirIn = [0.0, -0.707, -0.707] #[m]
        qMagIn = 10.0 #[W/m^2]
        # qMagIn = [10.0, 10.0] #[W/m^2]
#
        self.box = Solid.MeshSolid(stlPath, stpPath)
        #self.box = Solid.MeshSolid(stlPath, stpPath) #normally, use this one!

        self.fwd = ForwardModel.ForwardModel_MeshHF(self.box, qMagIn, qDirIn, hfMode='uniform') 
        # self.fwd = ForwardModel.ForwardModel_MeshHF(self.box, qMagIn, qDirIn, hfMode='uniform_multiple') 
        # self.fwd = ForwardModel.ForwardModel_MeshHF(self.box, qMagIn, qDirIn, hfMode='exponnorm') 
        self.opt = OptModel.OptModel_MeshHF()

        return
    

    def makeDirectories(self, runID, coefficientsList):
        c0 = coefficientsList[0]
        c1 = coefficientsList[1]
        c2 = coefficientsList[2]
        c3 = coefficientsList[3]
        # c4 = coefficientsList[4]
        runName = runID + f'_c0_{c0:.2f}_c1_{c1:.2f}_c2_{c2:.2f}_c3_{c3:.2f}' #f'_c0_{c0:.2f}_c1_{c1:.2f}_c2_{c2:.2f}_c3_{c3:.2f}'
        runName = runName.replace(".", "-")

        directoryName = f"{runName}" #running this within docker container means can't save to external without bindmount aaa

        os.makedirs(directoryName)
        # os.makedirs(f"{directoryName}/images")

        print(f"Made directory: {directoryName}")

        #where we are saving generated vtk's
        # directoryName = f"{directoryName}/vtks"
        # os.makedirs(directoryName)

        return directoryName

    def makeSweepCSV(self, c0List, c1List, c2List, c3List, maxhfList, runID):
        data = {
            'c0': c0List,
            'c1': c1List,
            'c2': c2List,
            'c3': c3List,
            'maxhf': maxhfList
        }
        df = pd.DataFrame(data)
        df.to_csv(f"{runID}_results.csv", index=None, header=True)

        return 
    
    # def findCornerFaces(self, trimeshSolid):
    #     # map each vertex to all the faces that contain it
    #     vertex_to_face_map = {}
    #     all_vertices = trimeshSolid.vertices
    #     for face_index, face in enumerate(trimeshSolid.faces):
    #         for vertex in face:
    #             if vertex in vertex_to_face_map:
    #                 vertex_to_face_map[vertex].append(face_index)
    #             else:
    #                 vertex_to_face_map[vertex] = [face_index]

    #     # identify corner faces
    #     corner_faces = set()
    #     # for vertex, faces in vertex_to_face_map.items():
    #     #     # print(vertex)
    #     #     # print(all_vertices[vertex])
    #     #     # if len(faces) !=6 and all_vertices[vertex][1] <= 10.0 and ((all_vertices[vertex][2] == 0.0) or (all_vertices[vertex][2] == 10.0)):  # a vertex is at a corner if it's not part of 6 faces, and we only want corners not the top edge 
    #     #     if len(faces) !=6 and ((all_vertices[vertex][1] == 10.0) or (all_vertices[vertex][1] == 0.0)) and ((all_vertices[vertex][2] == 0.0) or (all_vertices[vertex][2] == 10.0)):  # a vertex is at a corner if it's not part of 6 faces, and we only want corners not the top edge 
    #     #     #if len(faces) != 6:  # a vertex is at a corner if it's not part of 6 faces
    #     #         #should maybe add a condition for vertex y coordinates?
    #     #         #since this shouldn't stop us from moving elements on the top sharp edge - but we only want to 
    #     #         #exclude vertices on y=10.0, not ymax? 
    #     #         corner_faces.update(faces)
    #     return corner_faces

    def runOptimization(self):

        self.box.processSolid()
        trimeshSolid = self.box.trimeshSolid

        #find which faces are on corners - don't want to take those into account for difference in normals? 
        #reasoning is more than 1 adjacent face, and more than 1 normal to compare as a result? 
        #maybe this is not the way to go? 
        # corner_faces = self.findCornerFaces(self.box.trimeshSolid)
        adjacency = trimeshSolid.face_adjacency
        # filtered_adjacency = np.array([pair for pair in adjacency if pair[0] not in corner_faces and pair[1] not in corner_faces])
        
        def calculateNormalsDiff(trimeshSolid):
            """
            calculate differences between normal vectors of adjacent faces - ideally this should also be minimized?  
            also find max diff btwn adjacent normals
            """
            #"(n, 2) list of face indices. Each pair of faces in the list shares an edge, making them adjacent"
            # adjacency = trimeshSolid.face_adjacency 

            normals = trimeshSolid.face_normals
            normalsDiff = normals[adjacency[:, 0]] - normals[adjacency[:, 1]]
            normalsDiffMagnitude = np.linalg.norm(normalsDiff, axis=1)

            normals_0 = normals[adjacency[:, 0]]
            normals_1 = normals[adjacency[:, 1]]
            dot_product = np.einsum('ij,ij->i', normals_0, normals_1)
            clipped_dot_product = np.clip(dot_product, -1.0, 1.0)
            allAnglesBetweenNormals = np.arccos(clipped_dot_product)
            maxAngleBetweenNormals = np.max(allAnglesBetweenNormals)
            return normalsDiffMagnitude, maxAngleBetweenNormals
            # maxAngleDiff = np.pi - maxAngleBetweenNormals
            # return normalsDiffMagnitude, maxAngleDiff
            


        # def facesToKeep(mesh_center_xvals, mesh_center_yvals, mesh_center_zvals):

        #     # return np.where(
        #     #     #this should be for faces to NOT use
        #     #     #don't want to move anything below the slanted faces
        #     #     #don't want to move sides - keep those planar, and those are @ z=0 and z=10
        #     #     (mesh_center_yvals <= 10.0) #| #or should we keep this for sphere runs? 
        #     #     # (mesh_center_zvals == 0.0) |
        #     #     # (mesh_center_zvals == 10.0)
        #     # )[0]
        #     return [] #should we keep this for sphere runs? 

        def facesToKeep(trimeshSolid):
            #this should be for faces to NOT use (so they shouldn't move)
            return []
            # normals = trimeshSolid.face_normals
            # mesh_centers = trimeshSolid.triangles_center
            # mesh_center_yvals = mesh_centers[:, 1]
            # mesh_center_xvals = mesh_centers[:, 0]
            # mesh_center_zvals = mesh_centers[:, 2]
            # return np.where(
            #     #find normals in unit directions (+/- x, y, z)
            #     # (normals[:, 0] == 1.0) | (normals[:, 0] == -1.0) |    
            #     # (normals[:, 1] == 1.0) | (normals[:, 1] == -1.0) |
            #     # (normals[:, 2] == 1.0) | (normals[:, 2] == -1.0) 
                
            #     mesh_center_yvals <= 10.0
            # )[0]

        
        q_mesh_initial = self.fwd.calculateAllHF(trimeshSolid)
        maxHF_initial = self.fwd.filteredCalculateMaxHF(q_mesh_initial, unconstrainedFaces = [])
        sumHF_initial = self.fwd.calculateHFMeshSum(q_mesh_initial)
        # normalsDiff_initial, normalRefDotProducts_initial = calculateNormalsDiff(trimeshSolid)
        normalsDiff_initial, maxNormalsDiff_initial = calculateNormalsDiff(trimeshSolid)
        normalsPenalty_initial = np.sum(normalsDiff_initial)
        energy_initial = self.fwd.calculateIntegratedEnergy(q_mesh_initial, trimeshSolid)

        print(f"Initial max HF: {maxHF_initial}")
        print(f"Initial sum HF: {sumHF_initial}")
        print(f"Initial normals penalty: {normalsPenalty_initial}")
        print(f"Initial max angle between normals: {maxNormalsDiff_initial}")
        print(f"Initial energy: {energy_initial}")

        def objectiveFunction(trimeshSolid, coefficientsList, unconstrainedFaces):

            t0 = time.time()

            c0 = coefficientsList[0] #max heat flux term
            c1 = coefficientsList[1] #sum heat flux, unweighted, term
            c2 = coefficientsList[2] #normals diff term
            c3 = coefficientsList[3] #max normals term now 

            c4 = coefficientsList[4] #heat flux diff term

            q_mesh_all = self.fwd.calculateAllHF(trimeshSolid)

            # print(f"Time elapsed for q_mesh_all: {time.time() - t0}")

            unconstrainedFaces = [] #only keep this for cases with hf only on top face? 

            #for below: normalize all terms by dividing by initial values from pre-modification cube

            maxHFTerm = 0#c0 * (self.fwd.filteredCalculateMaxHF(q_mesh_all, unconstrainedFaces) / maxHF_initial)
            sumHFTerm = 0#c1 * (self.fwd.calculateHFMeshSum(q_mesh_all) / sumHF_initial) #self.fwd.calculateHFMeshSum(trimeshSolid)

            # normalsDiff, normalRefDotProducts = calculateNormalsDiff(trimeshSolid)
            normalsDiff, maxNormalsDiff = calculateNormalsDiff(trimeshSolid)
            normalsDiffSum = np.sum(normalsDiff)
            normalsPenalty = c2 * (normalsDiffSum / normalsPenalty_initial)
            maxNormalsTerm = c3 * maxNormalsDiff #(maxNormalsDiff / maxNormalsDiff_initial)

            # print(f"Normals diff sum: {normalsDiffSum}, max angle between normals: {maxNormalsDiff}")
            # input()

            # energyTerm = c3 * self.fwd.calculateIntegratedEnergy(q_mesh_all) #self.fwd.calculateIntegratedEnergy(trimeshSolid)
            # energyTerm = 0#c3 * (self.fwd.calculateIntegratedEnergy(q_mesh_all, trimeshSolid) / energy_initial)

            # print(f"Terms: {maxHFTerm}, {sumHFTerm}, {normalsPenalty}, {energyTerm}")
            # print(f"Terms divided by constants: {maxHFTerm/c0}, {sumHFTerm/c1}, {normalsPenalty/c2}, {energyTerm/c3}")

            # print(f"Time elapsed for whole objective calc: {time.time() - t0}")
            # input()

            return [maxHFTerm + sumHFTerm + normalsPenalty + maxNormalsTerm, normalsPenalty, maxNormalsTerm]
            # return maxHFTerm + sumHFTerm + normalsPenalty + energyTerm
            #return [maxHFTerm + sumHFTerm + normalsPenalty + energyTerm + hfDiffTerm, maxHFTerm, sumHFTerm, normalsPenalty, energyTerm] 
            #objectiveFunction value: [0]       
            

        def sweep_coefficients_and_record_output(coefficients_list, idx_to_vary, sweep_values):
            """
            run through possible values of one coefficient while keeping the rest fixed
            do optimization and save results and attempted combinations to csv
            """

            c0_runvals = []
            c1_runvals = []
            c2_runvals = []
            c3_runvals = []
            maxhf_vals = []

            for val in sweep_values:
                my_trimeshSolid = trimeshSolid.copy()
                coefficients_list[idx_to_vary] = val
                directoryName = self.makeDirectories(f"spheretest11_{idx_to_vary}", coefficients_list)
                #meshHFOpt(self, hfObjectiveFcn, constraint, updateHFProfile, calcHFAllMesh, calcMaxHF, calcEnergy, meshObj, coefficientsList, threshold, delta, id):
                maxHF = self.opt.meshHFOpt(
                    objectiveFunction,  
                    facesToKeep,
                    self.fwd.makeHFProfile,
                    self.fwd.calculateAllHF,
                    self.fwd.filteredCalculateMaxHF, #self.fwd.calculateMaxHF,
                    self.fwd.calculateIntegratedEnergy,
                    my_trimeshSolid, 
                    coefficients_list,
                    threshold=0.00000001, 
                    delta=0.01, 
                    id=directoryName
                )[0]
                c0_runvals.append(coefficients_list[0])
                c1_runvals.append(coefficients_list[1])
                c2_runvals.append(coefficients_list[2])
                c3_runvals.append(coefficients_list[3])
                maxhf_vals.append(maxHF)
            
            self.makeSweepCSV(c0_runvals, c1_runvals, c2_runvals, c3_runvals, maxhf_vals, f"{directoryName}/{idx_to_vary}")
            return 
        
        def big_sweep_coefficients_and_record_output(sweep_values_c0, sweep_values_c1, sweep_values_c2, sweep_values_c3, id):
            """
            run through possible values of one coefficient while keeping the rest fixed
            do optimization and save results and attempted combinations to csv
            """

            c0_runvals = []
            c1_runvals = []
            c2_runvals = []
            c3_runvals = []
            maxhf_vals = []

            for c0 in sweep_values_c0:
                for c1 in sweep_values_c1: 
                    for c2 in sweep_values_c2:
                        for c3 in sweep_values_c3: 
                            my_trimeshSolid = trimeshSolid.copy()
                            # coefficients_list[0] = c0
                            coefficients_list = [c0, c1, c2, c3, 0.0]
                            directoryName = self.makeDirectories(f"run{id}", coefficients_list)
                            #meshHFOpt(self, hfObjectiveFcn, constraint, updateHFProfile, calcHFAllMesh, calcMaxHF, calcEnergy, meshObj, coefficientsList, threshold, delta, id):
                            maxHF = self.opt.meshHFOpt(
                                objectiveFunction,  
                                facesToKeep,
                                self.fwd.makeHFProfile,
                                self.fwd.calculateAllHF,
                                self.fwd.filteredCalculateMaxHF, #self.fwd.calculateMaxHF,
                                self.fwd.calculateIntegratedEnergy,
                                my_trimeshSolid, 
                                coefficients_list,
                                threshold=0.00000001, 
                                delta=0.01, 
                                id=directoryName
                            )[0]
                            c0_runvals.append(coefficients_list[0])
                            c1_runvals.append(coefficients_list[1])
                            c2_runvals.append(coefficients_list[2])
                            c3_runvals.append(coefficients_list[3])
                            maxhf_vals.append(maxHF)
            
            self.makeSweepCSV(c0_runvals, c1_runvals, c2_runvals, c3_runvals, maxhf_vals, f"sweep3_{id}")
            return 

        ## For bulk variable sweep testing
        # sweep_values_c0 = [0, 50, 100, 150]
        # sweep_values_c1 = [0, 50, 100, 150]
        # sweep_values_c2 = [0, 50, 100, 150]
        # sweep_values_c3 = [0, 50, 100, 150]
        # big_sweep_coefficients_and_record_output(sweep_values_c0, sweep_values_c1, sweep_values_c2, sweep_values_c3, id=0)

        # sweep_values_c0 = [100, 200, 300, 400]
        # sweep_values_c1 = [100, 200, 300, 400]
        # sweep_values_c2 = [100, 200, 300, 400]
        # sweep_values_c3 = [100, 200, 300, 400]
        # big_sweep_coefficients_and_record_output(sweep_values_c0, sweep_values_c1, sweep_values_c2, sweep_values_c3, id=1)

        # sweep_values_c0 = [200, 350, 500]
        # sweep_values_c1 = [350, 400, 450]
        # sweep_values_c2 = [500, 1000, 1500]
        # sweep_values_c3 = [500, 2000, 4000]
        # big_sweep_coefficients_and_record_output(sweep_values_c0, sweep_values_c1, sweep_values_c2, sweep_values_c3, id=3)


        ## For variable sweep testing
        #coefficients_list = [21.16, 0.53, 14.0, 4.55, 0.0]
        # coefficients_list = [0.0, 100.0, 0.0, 0.0, 0.0]
        coefficients_list = [0, 0, 5000, 0, 0]
        # sweep_c0 = [0.0, 10.0, 20.0, 30.0]
        # sweep_c1 = [0.0, 0.53, 1.06, 1.59]
        # sweep_c2 = [0.0, 14.0, 28.0, 42.0]
        # sweep_c3 = [0.0, 2.275, 4.55, 6.825]

        # sweep_c0 = [1.0, 5.0, 10.0, 15.0]
        # sweep_c1 = [50, 100, 200, 500]
        # sweep_c2 = [50, 100, 200, 500]
        # sweep_c3 = [1.0, 5.0, 10.0, 15.0]

        # sweep_c3 = [300, 500, 1000, 2000]
        # sweep_c2 = [700, 1000, 2000, 5000]

        # sweep_c1 = [1.0, 10.0, 50.0, 100.0]
        # sweep_c1 = [50.0, 100.0, 300.0, 500.0]
        # sweep_c2 = [10.0, 50.0, 100.0, 300.0, 500.0, 1000.0, 1500.0, 2000.0]

        # sweep_c0 = [50.0, 100.0, 300.0, 500.0]

        # sweep_coefficients_and_record_output(coefficients_list, 0, sweep_c0)

        # sweep_coefficients_and_record_output(coefficients_list, 1, sweep_c1)

        # sweep_c2 = [4000, 5000, 6000]
        # sweep_c2 = [5000]
        # sweep_c2 = [100]
        # sweep_c2 = [3000, 3500, 4000, 4500, 5500, 6000]

        # sweep_coefficients_and_record_output(coefficients_list, 2, sweep_c2)

        # sweep_c3 = [1000, 5000, 10000]
        sweep_c3 = [250]
        sweep_coefficients_and_record_output(coefficients_list, 3, sweep_c3)

        # #more runs 
        # my_trimeshSolid = trimeshSolid.copy()
        # coefficientsList = [0.0, 100, 0.0, 0.0, 0.0]
        # # # coefficientsList = [21.16, 0.53, 8.95, c3]
        # # # coefficientsList = [0, 0, 0, c3]
        # # # coefficientsList = [30.0, 1.0, 14.0, 4.55]
        # # # coefficientsList = [21.16, 0.53, 14.0, 4.55]
        # # coefficientsList = [21.16, 0.53, 14.0, 4.55, 0.0]
        # # # coefficientsList = [1.0, 1.0, 1.0, 1.0, 1.0]
        # # # directoryName = self.makeDirectories(f"sweep_c4_{self.fwd.q_dir[0]}_{self.fwd.q_dir[1]}_{self.fwd.q_dir[2]}", coefficientsList)
        # directoryName = self.makeDirectories(f"2dir_test_5", coefficientsList)
        # # #meshHFOpt(self, hfObjectiveFcn, constraint, updateHFProfile, calcHFAllMesh, calcMaxHF, calcEnergy, meshObj, coefficientsList, threshold, delta, id):
        # self.opt.meshHFOpt(
        #         objectiveFunction,
        #         facesToKeep,  
        #         self.fwd.makeHFProfile,
        #         self.fwd.calculateAllHF,
        #         self.fwd.filteredCalculateMaxHF, #self.fwd.calculateMaxHF,
        #         self.fwd.calculateIntegratedEnergy,
        #         my_trimeshSolid, 
        #         coefficientsList,
        #         threshold=0.000001, 
        #         delta=0.01, 
        #         id=directoryName
        # )
        
        return 


if __name__ == '__main__':

    t0 = time.time()

    setup = RunSetup_MeshHF()

    """
    Run cube mesh operation test
    """
    setup.runOptimization()

    print(f"Time elapsed: {time.time() - t0}")





        
