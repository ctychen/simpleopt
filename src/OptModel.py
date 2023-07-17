import FreeCAD
import Part
import Mesh
import MeshPart
import numpy as np

import time

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import plotly.express as px

import pandas as pd
import multiprocessing

import vtk
from vtk import vtkPolyData, vtkPoints, vtkCellArray, vtkDoubleArray, vtkPolyDataWriter, vtkTriangle

from vtk.util.numpy_support import numpy_to_vtk
from scipy.optimize import curve_fit
from vtk.util import numpy_support

class OptModel_MeshHF: 
    """
    Mesh element manipulation to optimize heat flux - model for optimization
    """

    def __init__(self):
        """
        TODO: figure out what properties we actually need for optmodel, and what can just be fcn inputs
        """
        self.Ncores = multiprocessing.cpu_count() - 2 #reserve 2 cores for overhead
        #in case we run on single core machine
        if self.Ncores <= 0:
            self.Ncores = 1
        return
    
    # def gradientDescentHF(self, tri_mesh, objectiveFunction, allmeshelementsHF, facesToKeep, facesToMove, coefficientsList, delta, filedir, count):
    #     """
    #     gradient descent implementation for heat flux minimization
    #     takes in trimesh object and sorts elements by HF to deal with worst elements first
    #     calc gradient for each element by moving vertices a small amount and finding change in objective function
    #     move each vertex based on gradient * delta when all gradients calculated
    #     """ 

    #     use_set = set(np.where(allmeshelementsHF >= -10.0)[0]) #changed for 3sphere test
    #     gradient = np.zeros_like(tri_mesh.vertices)

    #     for idx in use_set:
    #         face = tri_mesh.faces[idx]
    #         for vertexIdx in face:  #vertexIdx is VERTEX INDICES
    #             obj_beforeMoving = objectiveFunction(tri_mesh, coefficientsList, facesToMove)[0]
    #             for j in range(3): #for every dimension - move the vertex a bit and calculate the change in objectiveFunction
    #                 tri_mesh.vertices[vertexIdx, j] += delta
    #                 obj_afterMoving = objectiveFunction(tri_mesh, coefficientsList, facesToMove)[0]
    #                 tri_mesh.vertices[vertexIdx, j] -= delta
    #                 gradient[vertexIdx, j] = (obj_afterMoving - obj_beforeMoving) / (2 * delta)
    #             #basically - move each vertex and update it
    #             tri_mesh.vertices[vertexIdx, 0] -= (delta * gradient[vertexIdx, 0])
    #             tri_mesh.vertices[vertexIdx, 1] -= (delta * gradient[vertexIdx, 1])
    #             tri_mesh.vertices[vertexIdx, 2] -= (delta * gradient[vertexIdx, 2])

    #     return tri_mesh


    def gradientDescentHF(self, tri_mesh, objectiveFunction, allmeshelementsHF, facesToKeep, facesToMove, coefficientsList, delta, filedir, count):
        f = tri_mesh.faces
        # v = tri_mesh.vertices
        v = tri_mesh.vertices.copy()
        gradient = np.zeros_like(tri_mesh.vertices)
        flattened_vertex_indices= f.flatten()
        reduced_vertices_indices = np.unique(flattened_vertex_indices) #and also apply any other 
        print(f"Reduced vertices: {reduced_vertices_indices}, length: {len(reduced_vertices_indices)}")
        # map reduced vertices back to original vertices
        map_from_reduced_vertices = np.array([np.where(tri_mesh.vertices == vertex)[0][0] for vertex in reduced_vertices_indices])
        print(f"Map from reduced vertices: {map_from_reduced_vertices}, length: {len(map_from_reduced_vertices)}")  

        for vertex_idx in reduced_vertices_indices:
            obj_beforeMoving = objectiveFunction(tri_mesh, coefficientsList, facesToMove)[0]
            for j in range(3):
                v[vertex_idx, j] += delta
                obj_afterMoving = objectiveFunction(tri_mesh, coefficientsList, facesToMove)[0]
                v[vertex_idx, j] -= delta   
                gradient[vertex_idx, j] = (obj_afterMoving - obj_beforeMoving) / (2 * delta)
  
            # v[vertex_idx, 0] -= (delta * gradient[vertex_idx, 0])
            # v[vertex_idx, 1] -= (delta * gradient[vertex_idx, 1])
            # v[vertex_idx, 2] -= (delta * gradient[vertex_idx, 2])

        #map back to original vertices - is this the right way to do this? 
        tri_mesh.vertices = v[map_from_reduced_vertices]

        return tri_mesh
    

    # def gradientDescentHF(self, tri_mesh, objectiveFunction, allmeshelementsHF, facesToKeep, facesToMove, coefficientsList, delta, filedir, count):
    #     """
    #     gradient descent implementation for heat flux minimization
    #     takes in trimesh object and sorts elements by HF to deal with worst elements first
    #     calc gradient for each element by moving vertices a small amount and finding change in objective function
    #     move each vertex based on gradient * delta when all gradients calculated
    #     """ 

    #     use_set = set(np.where(allmeshelementsHF >= -10.0)[0]) #changed for 3sphere test
    #     gradient = np.zeros_like(tri_mesh.vertices)

    #     def parallelGradientCalc(j):
    #         """
    #         calculate gradient for each vertex in parallel
    #         """
    #         tri_mesh.vertices[vertexIdx, j] += delta
    #         obj_afterMoving = objectiveFunction(tri_mesh, coefficientsList, facesToMove)[0]
    #         tri_mesh.vertices[vertexIdx, j] -= delta
    #         gradientValue = (obj_afterMoving - obj_beforeMoving) / (2 * delta)
    #         return gradientValue

    #     for idx in use_set:
    #         face = tri_mesh.faces[idx]
    #         for vertexIdx in face:  #vertexIdx is VERTEX INDICES
    #             obj_beforeMoving = objectiveFunction(tri_mesh, coefficientsList, facesToMove)[0]
    #             # for j in range(3): #for every dimension - move the vertex a bit and calculate the change in objectiveFunction
    #             #     tri_mesh.vertices[vertexIdx, j] += delta
    #             #     obj_afterMoving = objectiveFunction(tri_mesh, coefficientsList, facesToMove)[0]
    #             #     tri_mesh.vertices[vertexIdx, j] -= delta
    #             #     gradient[vertexIdx, j] = (obj_afterMoving - obj_beforeMoving) / (2 * delta)
    #             #basically - move each vertex and update it
    #             pool = multiprocessing.Pool(self.Ncores)
    #             gradient[vertexIdx, :] = np.asarray(pool.map(parallelGradientCalc, np.arange(3)))  
    #             print(f"gradient: {gradient}")
    #             pool.close()
    #             pool.join() 
    #             del pool

    #             tri_mesh.vertices[vertexIdx, 0] -= (delta * gradient[vertexIdx, 0])
    #             tri_mesh.vertices[vertexIdx, 1] -= (delta * gradient[vertexIdx, 1])
    #             tri_mesh.vertices[vertexIdx, 2] -= (delta * gradient[vertexIdx, 2])

    #     return tri_mesh


    # def meshHFOpt(self, hfObjectiveFcn, calcHFAllMesh, calcMaxHF, calcEnergy, meshObj, coefficientsList, threshold, delta, id):
    def meshHFOpt(self, hfObjectiveFcn, constraint, updateHFProfile, calcHFAllMesh, calcMaxHF, calcEnergy, meshObj, coefficientsList, threshold, delta, id):
    # def meshHFOpt(self, hfFunction, hfObjectiveFcn, meshObj, threshold, step, id):
        """
        runs optimization process until objective fcn value reaches stopping condition @ minimum
        modifies the mesh based on gradient by applying changeMeshFcn accordingly

        can change changeMeshFcn, hfObjectiveFcn to use different functions
        if we want a different manipulation, or add more stuff to the functions

        can also set constraint to be whatever conditions should be true for the faces we can manipulate. 
        basically, if the constraint is true, we can move the face, otherwise we won't do anything to it
        """
        #TODO: add constraints somehow - take in list of criteria? eg. don't move face if x=0 or y=0 or x=10 or y=10?

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
        objFcn = hfObjectiveFcn(trimeshSolid, coefficientsList, facesToMove)
        # all_objective_function_values = [hfObjectiveFcn(trimeshSolid, coefficientsList, facesToMove)]
        all_objective_function_values = [objFcn[0]]
        all_sum_normals_diff = [objFcn[1]]
        all_max_normals_diff = [objFcn[2]]
        hf_all_mesh = calcHFAllMesh(trimeshSolid)
        max_hf_each_run = [calcMaxHF(hf_all_mesh, facesToMove)]

        # make VTK to display HF on surface
        self.plotHFVTK(calcHFAllMesh(trimeshSolid), trimeshSolid, f"{id}", count=4242)

        prev_objVal = 2000
        curr_objVal = 0

        # t0 = time.time()

        #faces to NOT move
        facesToKeep = indicesToNotMove

        while abs(prev_objVal - curr_objVal) > threshold and count < 500:

            hf_all_mesh = calcHFAllMesh(trimeshSolid)

            #calc the gradient
            trimeshSolid = self.gradientDescentHF(trimeshSolid, hfObjectiveFcn, hf_all_mesh, facesToKeep, facesToMove, coefficientsList, delta, f"test{id}", count)
            
            #recalculate the hf profile on the surface - don't need this for spheretests
            # updateHFProfile(trimeshSolid) 

            # new_objVal = hfObjectiveFcn(trimeshSolid, coefficientsList, facesToMove)
            new_objVal_all = hfObjectiveFcn(trimeshSolid, coefficientsList, facesToMove)
            new_objVal = new_objVal_all[0]
            new_sum_normals_diff = new_objVal_all[1]
            new_max_normals_diff = new_objVal_all[2]
            all_objective_function_values.append(new_objVal)
            all_sum_normals_diff.append(new_sum_normals_diff)
            all_max_normals_diff.append(new_max_normals_diff)

            prev_objVal = curr_objVal
            curr_objVal = new_objVal

            new_max_hf = calcMaxHF(hf_all_mesh, facesToMove)
            max_hf_each_run.append(new_max_hf)

            ##this is for plotting surface fit onto mesh
            # if count % 5 == 0:
            #     self.makePolyFitSurface(trimeshSolid.vertices[constrainedVIdx], f"{id}", count)

            # print(f"New objective function value: {new_objVal}")

            if count and count % 10 == 0: #count % 5 == 0: 

                # x_count = np.linspace(0, len(max_hf_each_run), len(max_hf_each_run))
                # fig = px.scatter(x = x_count, y = max_hf_each_run)
                # fig.update_xaxes(title_text='Iterations')
                # fig.update_yaxes(title_text=f'{calcMaxHF.__name__}')
                # fig.show()            
                # output_file = f"{id}/max_hf_run_{count}.html"
                # pio.write_html(fig, output_file)

                # x_count = np.linspace(0, len(all_objective_function_values), len(all_objective_function_values))
                # fig = px.scatter(x = x_count, y = all_objective_function_values)
                # fig.update_xaxes(title_text='Iterations')
                # fig.update_yaxes(title_text='Objective function - sum HF over elements')
                # fig.show()            
                # output_file = f"{id}/objective_run_{count}.html"
                # pio.write_html(fig, output_file)

                # #make VTK to display HF on surface
                self.plotHFVTK(calcHFAllMesh(trimeshSolid), trimeshSolid, f"{id}", count)

            # if count % 10 == 0: 
                # self.plotHFVTK(calcHFAllMesh(trimeshSolid), trimeshSolid, f"{id}", count)
                # self.plotHFVTK(calcHFAllMesh(trimeshSolid), trimeshSolid, f"test{id}")

            count += 1
        
        # self.plotRun(all_objective_function_values, max_hf_each_run, f"{id}")
        self.plotHFVTK(calcHFAllMesh(trimeshSolid), trimeshSolid, f"{id}", count)
        self.plotMaxNormalsDiff(all_max_normals_diff, f"{id}")
        self.plotNormalsDiff(all_sum_normals_diff, f"{id}")
        self.plotObjectiveFunction(all_objective_function_values, f"{id}")
        self.plotMaxHF(max_hf_each_run, f"{id}")    
        # finalMaxHF = np.min(max_hf_each_run)
        # print(f"Finished run, maxHF is {finalMaxHF}")
        print(f"Finished run")
        
        return trimeshSolid
        # return finalMaxHF, trimeshSolid

        #when process is done, the mesh should have been modified - so return it 
        # return trimeshSolid

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