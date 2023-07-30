import FreeCAD
import Part
import Mesh
import MeshPart
import numpy as np

import time
import copy

import Solid
import ForwardModel
import ObjectiveFunctionTools
import trimesh

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import plotly.express as px

import pandas as pd
import multiprocessing
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor

import vtk
from vtk import vtkPolyData, vtkPoints, vtkCellArray, vtkDoubleArray, vtkPolyDataWriter, vtkTriangle

from vtk.util.numpy_support import numpy_to_vtk
from scipy.optimize import curve_fit
from vtk.util import numpy_support

### GRADIENT DESCENT AND OPTIMIZATION###

objfcnTools = ObjectiveFunctionTools.ObjectiveFunctionTools()

class OptModel_MeshHF: 
    """
    Mesh element manipulation to optimize heat flux - model for optimization
    """

    def __init__(self):
        """
        TODO: figure out what properties we actually need for optmodel, and what can just be fcn inputs
        """
        self.Ncores = 2 #multiprocessing.cpu_count() - 1
        #in case we run on single core machine
        if self.Ncores <= 0:
            self.Ncores = 1

        print(f"number cores being used: {self.Ncores}")
        return

    #gradientDescentHF(self, tri_mesh, allmeshelementsHF, face_adjacency, face_adjacency_edges, initialParams, facesToKeep, facesToMove, coefficientsList, delta, filedir, count):
    def gradientDescentHF(self, tri_mesh, allmeshelementsHF, facesToMove, delta):

        use_set = set(np.where(allmeshelementsHF >= -10.0)[0]) #changed for 3sphere test
        vertices = tri_mesh.vertices
        gradient = np.zeros_like(vertices)

        all_faces = tri_mesh.faces

        useFaces = all_faces[list(use_set)]
        objfcnTools.setFacesToMove(facesToMove)
        flattenedVtx = useFaces.flatten() 
        uniqueVtx = np.unique(flattenedVtx)
        numVtx = len(vertices)

        currentVerticesGrid = np.array([np.tile(vertices[np.newaxis, :], (numVtx, 1, 1)) for _ in range(3)])
        numDim = len(currentVerticesGrid)

        currentObjFcnVal = objfcnTools.vtxFacesObjectiveFunctionCalc(currentVerticesGrid[0][0])
        currentObjectiveFcnValues = np.full((numVtx, 3), currentObjFcnVal)

        newVerticesGrid = currentVerticesGrid.copy()
        delta = delta * (254.0 / numVtx)
        range_indices = np.arange(currentVerticesGrid.shape[1])
        newVerticesGrid[0, range_indices, range_indices, 0] += delta
        newVerticesGrid[1, range_indices, range_indices, 1] += delta
        newVerticesGrid[2, range_indices, range_indices, 2] += delta
        objfcnTools.setNewVerticesGrid(newVerticesGrid)

        numProcesses = self.Ncores

        try: 
            pool = multiprocessing.Pool(numProcesses)
            newObjectiveFcnValues = np.array(pool.starmap(
                objfcnTools.objectiveFunction, 
                [(vtx, dim) for vtx in range(numVtx) for dim in range(numDim)]
            )).reshape(numVtx, 3)
        finally:
            pool.close()
            pool.join()
            del pool

        gradient = (newObjectiveFcnValues - currentObjectiveFcnValues) / (2 * delta) 

        # print(f"calculated new gradient")

        #basically - move each vertex and update it
        tri_mesh.vertices[uniqueVtx, 0] -= (delta * gradient[uniqueVtx, 0])
        tri_mesh.vertices[uniqueVtx, 1] -= (delta * gradient[uniqueVtx, 1])
        tri_mesh.vertices[uniqueVtx, 2] -= (delta * gradient[uniqueVtx, 2])

        return tri_mesh

    # def meshHFOpt(self, hfObjectiveFcn, calcHFAllMesh, calcMaxHF, calcEnergy, meshObj, coefficientsList, threshold, delta, id):
    def meshHFOpt(self, constraint, updateHFProfile, meshObj, coefficientsList, threshold, delta, fwdModel, id):
    # def meshHFOpt(self, hfFunction, hfObjectiveFcn, meshObj, threshold, step, id):
        """
        runs optimization process until objective fcn value reaches stopping condition @ minimum
        modifies the mesh based on gradient by applying changeMeshFcn accordingly
        """

        #assuming input is already a trimesh, ie. processing the solid was done already
        trimeshSolid = meshObj

        count = 0
 
        # #commented out for spheretests since don't need this 
        mesh_centers = trimeshSolid.triangles_center
        mesh_center_yvals = mesh_centers[:, 1]
        mesh_center_xvals = mesh_centers[:, 0]
        mesh_center_zvals = mesh_centers[:, 2]

        # indicesToNotMove = set(constraint(mesh_center_xvals, mesh_center_yvals, mesh_center_zvals))
        indicesToNotMove = set(constraint(trimeshSolid))
        allIndices = set(range(len(mesh_center_xvals)))  # assuming all arrays are the same length
        facesToMove = allIndices - indicesToNotMove #faces to move

        print(f"Objective function with coefficients: {coefficientsList}")
        vertices = trimeshSolid.vertices
        # faces = trimeshSolid.faces
        face_adjacency = trimeshSolid.face_adjacency
        face_adjacency_edges = trimeshSolid.face_adjacency_edges
        initialVolume = trimeshSolid.volume
        print(f"Initial volume: {initialVolume}")
        initialParams = [initialVolume] 

        #objective function tools setup 
        objfcnTools.setParams(initialParams, coefficientsList)
        objfcnTools.setForwardModel(fwdModel)
        objfcnTools.setMeshAndGrids(trimeshSolid)

        objFcnVal = objfcnTools.vtxFacesObjectiveFunctionCalc(vertices)
        all_objective_function_values = [objFcnVal]

        hf_all_mesh = objfcnTools.calculateAllHF(fwdModel.hfMode, fwdModel.q_dir, fwdModel.q_mag, vertices) #
        # max_hf_each_run = [calcMaxHF(hf_all_mesh, facesToMove)]

        # make VTK to display HF on surface
        self.plotHFVTK(objfcnTools.calculateAllHF(fwdModel.hfMode, fwdModel.q_dir, fwdModel.q_mag, vertices), trimeshSolid, f"{id}", count=4242)

        prev_objVal = 2000
        curr_objVal = 0

        #faces to NOT move
        facesToKeep = indicesToNotMove

        while abs(prev_objVal - curr_objVal) > threshold and count < 2000:

            vertices = trimeshSolid.vertices

            hf_all_mesh = objfcnTools.calculateAllHF(fwdModel.hfMode, fwdModel.q_dir, fwdModel.q_mag, vertices)

            t1 = time.time()

            #calculate the gradient
            trimeshSolid = self.gradientDescentHF(trimeshSolid, hf_all_mesh, facesToMove, delta)

            print(f"{count}: gradient descent time: {time.time() - t1}")

            vertices = trimeshSolid.vertices

            #recalculate the hf profile on the surface - don't need this for spheretests
            # updateHFProfile(trimeshSolid) 

            new_objVal = objfcnTools.vtxFacesObjectiveFunctionCalc(vertices)
            all_objective_function_values.append(new_objVal)

            prev_objVal = curr_objVal
            curr_objVal = new_objVal

            if count and count % 20 == 0: #count % 5 == 0: 
                #self.plotHFVTK(calcHFAllMesh(trimeshSolid), trimeshSolid, f"{id}", count)
                self.plotHFVTK(objfcnTools.calculateAllHF(fwdModel.hfMode, fwdModel.q_dir, fwdModel.q_mag, vertices), trimeshSolid, f"{id}", count)

            count += 1

        vertices = trimeshSolid.vertices
        
        # self.plotRun(all_objective_function_values, max_hf_each_run, f"{id}")
        # self.plotHFVTK(calcHFAllMesh(trimeshSolid), trimeshSolid, f"{id}", count)
        self.plotHFVTK(objfcnTools.calculateAllHF(fwdModel.hfMode, fwdModel.q_dir, fwdModel.q_mag, vertices), trimeshSolid, f"{id}", count)
        # self.plotMaxNormalsDiff(all_max_normals_diff, f"{id}")
        # self.plotNormalsDiff(all_sum_normals_diff, f"{id}")
        self.plotObjectiveFunction(all_objective_function_values, f"{id}")
        # self.plotMaxHF(max_hf_each_run, f"{id}")  
          
        # finalMaxHF = np.min(max_hf_each_run)
        # print(f"Finished run, maxHF is {finalMaxHF}")
        print(f"Finished run")
        
        return trimeshSolid

    def plotMaxNormalsDiff(self, max_normals_diff_runs, directoryName):
        x_count = np.linspace(0, len(max_normals_diff_runs), len(max_normals_diff_runs))
        fig = px.scatter(x = x_count, y = max_normals_diff_runs)
        fig.update_xaxes(title_text='Iterations')
        fig.update_yaxes(title_text='Max normals difference')
        fig.show()            
        output_file = f"{directoryName}/max_normals_diff_each_run.html"
        pio.write_html(fig, output_file)

    def plotNormalsDiff(self, normals_diff_runs, directoryName):
        x_count = np.linspace(0, len(normals_diff_runs), len(normals_diff_runs))
        fig = px.scatter(x = x_count, y = normals_diff_runs)
        fig.update_xaxes(title_text='Iterations')
        fig.update_yaxes(title_text='Sum of normals difference')
        fig.show()            
        output_file = f"{directoryName}/sum_normals_diff_each_run.html"
        pio.write_html(fig, output_file)

    def plotObjectiveFunction(self, objectiveFunctionValues, directoryName):
        x_count = np.linspace(0, len(objectiveFunctionValues), len(objectiveFunctionValues))
        fig = px.scatter(x = x_count, y = objectiveFunctionValues)
        fig.update_xaxes(title_text='Iterations')
        fig.update_yaxes(title_text='Objective Function Values')
        fig.show()            
        output_file = f"{directoryName}/objective_each_run.html"
        pio.write_html(fig, output_file)

    def plotMaxHF(self, maxHFValues, directoryName):
        x_count = np.linspace(0, len(maxHFValues), len(maxHFValues))
        fig = px.scatter(x = x_count, y = maxHFValues)
        fig.update_xaxes(title_text='Iterations')
        fig.update_yaxes(title_text='Max HF Values')
        fig.show()            
        output_file = f"{directoryName}/max_hf_each_run.html"
        pio.write_html(fig, output_file)


    def plotHFVTK(self, hfValues, trimeshSolid, fileDir, count):
        """
        Make and export VTK for visualizing HF on mesh elements for each iteration
        """

        # Now, we'll use VTK to create a colored mesh
        points = vtkPoints()
        polys = vtkCellArray()
        heatFluxMagnitudes = vtkDoubleArray()

        # Add the vertices
        for vertex in trimeshSolid.vertices:
            points.InsertNextPoint(vertex)

        # Add the faces (triangles)
        for face in trimeshSolid.faces:
            triangle = vtkTriangle()
            triangle.GetPointIds().SetId(0, face[0])
            triangle.GetPointIds().SetId(1, face[1])
            triangle.GetPointIds().SetId(2, face[2])
            polys.InsertNextCell(triangle)

        # Add the colors
        for value in hfValues:
            heatFluxMagnitudes.InsertNextValue(value)

        # Create the mesh object
        polydata = vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetPolys(polys)
        polydata.GetCellData().SetScalars(heatFluxMagnitudes)

        # Write to a .vtk file
        writer = vtkPolyDataWriter()
        writer.SetFileName(f"{fileDir}/{count:05}.vtk")
        writer.SetInputData(polydata)
        writer.Write()

        return 
    
    def makePolyFitSurface(self, verticesToUse, fileDir, count):

        def poly_surface(coord, a, b, c, d, e, f):
            x, y = coord
            return a * x**2 + b * y**2 + c * x * y + d * x + e * y + f
        
        # Create a 2D coordinate grid
        x = verticesToUse[:, 0]
        z = verticesToUse[:, 2]
        coord = np.vstack((x, z))

        # y values are the heights for the surface
        y = verticesToUse[:, 1]

        # Fit the surface using SciPy's curve_fit function
        popt, pcov = curve_fit(poly_surface, coord, y)

        # Calculate distances from the points to the fitted surface
        y_fit = poly_surface(coord, *popt)
        distances = np.abs(y - y_fit)

        # Visualization using VTK

        # Create a structured grid
        grid = vtk.vtkStructuredGrid()

        # Generate a grid over the x and z values
        x_grid, z_grid = np.meshgrid(np.linspace(x.min(), x.max(), 100), np.linspace(z.min(), z.max(), 100))

        # Compute y values over the grid
        #curvefit returns values of params, so [a, b, c, d, e, f] coefficients
        y_fit_grid = poly_surface((x_grid.ravel(), z_grid.ravel()), *popt)

        # Combine x, y, z arrays to a single 3D array for the VTK grid
        points = np.column_stack([x_grid.ravel(), y_fit_grid, z_grid.ravel()])
        points_vtk = numpy_to_vtk(points, deep=1)

        # Define points and grid
        points = vtk.vtkPoints()
        points.SetData(points_vtk)
        grid.SetPoints(points)

        # Set dimensions of the grid
        grid.SetDimensions(x_grid.shape[0], x_grid.shape[1], 1)

        # Save grid to a .vts file
        writer = vtk.vtkXMLStructuredGridWriter()
        writer.SetFileName(f"{fileDir}/polyfit_{count:05}.vts")
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(grid)
        else:
            writer.SetInputData(grid)
        writer.Write()

        return 
    

    def plotRun(self, objective_function_values, max_hf_each_run, outputDir):
        """
        plot values of objective function, as well as max HF and sum of HF's, over iterations
        """
        x_count = np.linspace(0, len(objective_function_values), len(objective_function_values))
        fig = px.scatter(x = x_count, y = objective_function_values)
        fig.update_xaxes(title_text='Iterations')
        fig.update_yaxes(title_text='Objective function - sum HF over elements')
        fig.show()            
        output_file = f"{outputDir}/entire_run.html"
        pio.write_html(fig, output_file)

        x_count = np.linspace(0, len(max_hf_each_run), len(max_hf_each_run))
        fig = px.scatter(x = x_count, y = max_hf_each_run)
        fig.update_xaxes(title_text='Iterations')
        fig.update_yaxes(title_text='Max HF')
        fig.show()            
        output_file = f"{outputDir}/max_hf_each_run.html"
        pio.write_html(fig, output_file)

        return 



class OptModel_Template:
    """
    Template class for optimizer, change these functions to match problem we're dealing with 
    """
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