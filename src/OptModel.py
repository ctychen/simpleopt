import FreeCAD
import Part
import Mesh
import MeshPart
import numpy as np

import time
import copy

import Solid
import ForwardModel
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


# ###OBJECTIVE FUNCTIONS AND CONSTRAINTS###

# # def objective_for_vertex_dim(objectiveFunction, newVerticesGrid, vtx, dim, all_faces, coefficientsList, facesToMove):
# #     return objectiveFunction(trimesh.Trimesh(vertices=newVerticesGrid[dim][vtx], faces=all_faces), coefficientsList, facesToMove)[0]

# def objective_for_vertex_dim(objectiveFunction, newVerticesGrid, vtx, dim, all_faces, face_adjacency, face_adjacency_edges, initialParams, coefficientsList, facesToMove):
#     #face_adjacency, face_adjacency_edges, initialParams
#     return objectiveFunction(newVerticesGrid[dim][vtx], all_faces, face_adjacency, face_adjacency_edges, initialParams, coefficientsList, facesToMove)[0]

# # TODO - eventually when we need to call forward model functions - would need this 
# # but maybe better solution is to not have a FwdModel class/object and just have a bunch of functions in a file that we can call
# # def objective_for_vertex_dim(objectiveFunction, newVerticesGrid, vtx, dim, all_faces, coefficientsList, facesToMove, fwd):
# #     return objectiveFunction(trimesh.Trimesh(vertices=newVerticesGrid[dim][vtx], faces=all_faces), coefficientsList, facesToMove, fwd)

# def calculateNormalsDiff(trimeshSolid):
#     vertex_defects = trimesh.curvature.vertex_defects(trimeshSolid)
#     sumVertexDefects = np.sum(np.abs(vertex_defects))
#     maxVertexDefect = np.max(np.abs(vertex_defects))  
#     maxAngleBetweenNormals = 0 #TODO - will need to bring in these terms from impelmentation in RunModel
#     return sumVertexDefects, maxVertexDefect, maxAngleBetweenNormals

# def calculateAngles(pairs):
#     pairs = np.asanyarray(pairs, dtype=np.float64)
#     # do the dot product between vectors
#     dots = np.dot(pairs[:, 0] * pairs[:, 1], [1.0] * pairs[:, 0].shape[1])
#     # clip for floating point error
#     dots = np.clip(dots, -1.0, 1.0)
#     # do cos and remove arbitrary sign
#     angles = np.abs(np.arccos(dots))
#     return angles

# def calculateDotProducts(vecA, vecB):
#     #from trimesh - util.py
#     vecA = np.asanyarray(vecA)
#     return np.dot(vecA * vecB, [1.0] * vecA.shape[1])

# def calculateIntegralMeanCurvature(vertices, faces, face_adjacency, face_adjacency_edges):
#     #calc face normals
#     face_normals = Solid.calculateFaceNormals(vertices, faces)
#     #calc angles between adjacent faces
#     pairs = face_normals[face_adjacency]
#     angles = calculateAngles(pairs)
#     #calc integral mean curvature
#     edges_length = np.linalg.norm(np.subtract(
#         *vertices[face_adjacency_edges.T]), axis=1)
#     integralMeanCurvature = (angles * edges_length).sum() * 0.5
#     return integralMeanCurvature

# def calcUnitVectors(vectors):
#     norm = np.sqrt(np.dot(vectors, vectors))
#     unitVectors = vectors / norm
#     return unitVectors

# def calculateVertexDefects(vertices, faces, face_adjacency):
#     #calculate mesh face angles - so angles of triangles
#     triangles = np.asanyarray(vertices.view(np.ndarray)[faces], dtype=np.float64)
#     # get a unit vector for each edge of the triangle
#     u = calcUnitVectors(triangles[:, 1] - triangles[:, 0])
#     v = calcUnitVectors(triangles[:, 2] - triangles[:, 0])
#     w = calcUnitVectors(triangles[:, 2] - triangles[:, 1])

#     # run the cosine and per-row dot product
#     angles = np.zeros((len(triangles), 3), dtype=np.float64)
#     # clip to make sure we don't float error past 1.0
#     angles[:, 0] = np.arccos(np.clip(calculateDotProducts(u, v), -1, 1))
#     angles[:, 1] = np.arccos(np.clip(calculateDotProducts(-u, w), -1, 1))
#     # the third angle is just the remaining
#     angles[:, 2] = np.pi - angles[:, 0] - angles[:, 1]

#     # a triangle with any zero angles is degenerate
#     # so set all of the angles to zero in that case
#     angles[(angles < 1e-8).any(axis=1), :] = 0.0

#     angle_sum = np.array(angles.sum(axis=1)).flatten()

#     return (2*np.pi) - angle_sum


# #face_adjacency, faces, and face_adjacency_edges are from trimesh but all will only need to be accessed once at beginning, and are all np arrays
# def objectiveFunction(vertices, faces, face_adjacency, face_adjacency_edges, initialParams, coefficientsList, unconstrainedFaces):
#     #for initialParams - this is volume of original mesh, etc. - for normalization, etc. 
#     c0, c1, c2, c3, c4 = coefficientsList
#     # maxHFTerm = 0 #c0 * self.fwd.filteredCalculateMaxHF(q_mesh_all, unconstrainedFaces)    #try not dividing by initial value
#     # sumHFTerm = 0 #c1 * (self.fwd.calculateHFMeshSum(q_mesh_all) / numFaces) 
#     # sumVertexDefects, maxVertexDefects, maxAngleBetweenNormals = calculateNormalsDiff(trimeshSolid)  
#     imcTerm = c2 * calculateIntegralMeanCurvature(vertices, faces, face_adjacency, face_adjacency_edges)
#     # imcTerm = c2 * (Solid.calculateSurfaceArea(vertices, faces) / initialParams[0])
#     # normalsPenalty = c2 * sumVertexDefects
#     vertexDefectsTerm = 0#c2 * calculateVertexDefects(vertices, faces, face_adjacency)
#     maxNormalsTerm = 0#c3 * maxAngleBetweenNormals   
#     #c4 was originally a thing but i've given up
#     # return [vertexDefectsTerm + maxNormalsTerm, vertexDefectsTerm, maxNormalsTerm]
#     #return [maxHFTerm + sumHFTerm + normalsPenalty + maxNormalsTerm + maxVertexDefectsTerm, normalsPenalty, maxNormalsTerm, maxVertexDefectsTerm]
#     # return [normalsPenalty + maxNormalsTerm + maxVertexDefectsTerm, normalsPenalty, maxNormalsTerm]   
#     return [imcTerm + maxNormalsTerm + vertexDefectsTerm, imcTerm, maxNormalsTerm] 


### GRADIENT DESCENT AND OPTIMIZATION###

class OptModel_MeshHF: 
    """
    Mesh element manipulation to optimize heat flux - model for optimization
    """

    def __init__(self):
        """
        TODO: figure out what properties we actually need for optmodel, and what can just be fcn inputs
        """
        self.Ncores = multiprocessing.cpu_count()
        #in case we run on single core machine
        if self.Ncores <= 0:
            self.Ncores = 1
        return
    

    # def gradientDescentHF(self, tri_mesh, objectiveFunction, allmeshelementsHF, facesToKeep, facesToMove, coefficientsList, delta, filedir, count):
    # #     """
    # #     gradient descent implementation for heat flux minimization
    # #     takes in trimesh object and sorts elements by HF to deal with worst elements first
    # #     calc gradient for each element by moving vertices a small amount and finding change in objective function
    # #     move each vertex based on gradient * delta when all gradients calculated
    # #     """ 
    #     use_set = set(np.where(allmeshelementsHF >= -10.0)[0]) #changed for 3sphere test
    #     gradient = np.zeros_like(tri_mesh.vertices)
    #     useFaces = tri_mesh.faces[list(use_set)]
    #     flattenedVtx = useFaces.flatten() 
    #     uniqueVtx = np.unique(flattenedVtx)
    #     obj_beforeMoving = objectiveFunction(tri_mesh, coefficientsList, facesToMove)[0]
    #     for vertex in np.nditer(uniqueVtx):
    #         for j in range(3): 
    #             tri_mesh.vertices[vertex, j] += delta
    #             obj_afterMoving = objectiveFunction(tri_mesh, coefficientsList, facesToMove)[0]
    #             tri_mesh.vertices[vertex, j] -= delta
    #             gradient[vertex, j] = (obj_afterMoving - obj_beforeMoving) / (2 * delta)
    #     #basically - move each vertex and update it
    #     tri_mesh.vertices[uniqueVtx, 0] -= (delta * gradient[uniqueVtx, 0])
    #     tri_mesh.vertices[uniqueVtx, 1] -= (delta * gradient[uniqueVtx, 1])
    #     tri_mesh.vertices[uniqueVtx, 2] -= (delta * gradient[uniqueVtx, 2])

    #     return tri_mesh


### let's try loop-free approach... ###

    def gradientDescentHF(self, tri_mesh, objectiveFunction, allmeshelementsHF, face_adjacency, face_adjacency_edges, initialParams, facesToKeep, facesToMove, coefficientsList, delta, filedir, count):

        use_set = set(np.where(allmeshelementsHF >= -10.0)[0]) #changed for 3sphere test
        gradient = np.zeros_like(tri_mesh.vertices)

        all_faces = tri_mesh.faces

        useFaces = all_faces[list(use_set)]
        flattenedVtx = useFaces.flatten() 
        uniqueVtx = np.unique(flattenedVtx)
        numVtx = len(tri_mesh.vertices)

        currentVerticesGrid = np.array([np.tile(tri_mesh.vertices[np.newaxis, :], (numVtx, 1, 1)) for _ in range(3)])
        #currentObjFcnVal = objectiveFunction(trimesh.Trimesh(vertices=currentVerticesGrid[0][0], faces=all_faces), coefficientsList, facesToMove)[0]
        # currentObjFcnVal = objectiveFunction(currentVerticesGrid[0][0], all_faces, coefficientsList, facesToMove)[0]
        #objectiveFunction(vertices, faces, face_adjacency, face_adjacency_edges, coefficientsList, unconstrainedFaces):

        currentObjFcnVal = objectiveFunction(currentVerticesGrid[0][0], all_faces, face_adjacency, face_adjacency_edges, initialParams, coefficientsList, facesToMove)[0]
        currentObjectiveFcnValues = np.full((numVtx, 3), currentObjFcnVal)

        newVerticesGrid = currentVerticesGrid.copy()
        range_indices = np.arange(currentVerticesGrid.shape[1])
        newVerticesGrid[0, range_indices, range_indices, 0] += delta
        newVerticesGrid[1, range_indices, range_indices, 1] += delta
        newVerticesGrid[2, range_indices, range_indices, 2] += delta

        numProcesses = self.Ncores
        # someObject.newVerticesGrid = newVerticesGrid


        with Pool(numProcesses) as p:
            newObjectiveFcnValues = np.array(p.starmap(
                objective_for_vertex_dim, 
                [(objectiveFunction, newVerticesGrid, vtx, dim, all_faces, face_adjacency, face_adjacency_edges, initialParams, coefficientsList, facesToMove) for vtx in range(numVtx) for dim in range(len(newVerticesGrid))]
            )).reshape(numVtx, 3)

        # try:
        #     pool = multiprocessing.Pool(Ncores)
        #     #Moller-Trumbore algorithm
        #     if mode == 'MT':
        #         print("Using Moller-Trumbore intersection algorithm")
        #         log.info("Using Moller-Trumbore intersection algorithm")
        #         mask = np.asarray(pool.map(tools.intersectTestParallelMT, np.arange(N)))
        #     #Signed Volume algorithm
        #     else:
        #         print("Using signed volume intersection algorithm")
        #         log.info("Using signed volume intersection algorithm")
        #         mask = np.asarray(pool.map(tools.intersectTestParallel, np.arange(N)))
        # finally:
        #     pool.close()
        #     pool.join()
        #     del pool

        #def objective_for_vertex_dim(objectiveFunction, newVerticesGrid, vtx, dim, all_faces, face_adjacency, coefficientsList, facesToMove):
        #objectiveFunction(vertices, faces, face_adjacency, face_adjacency_edges, coefficientsList, unconstrainedFaces):

        # print(f"time to calculate NEW objfcn: {time.time() - t0}")

        gradient = (newObjectiveFcnValues - currentObjectiveFcnValues) / (2 * delta) 

        # print(f"gradient: {gradient}")
        delta = delta * (254.0 / numVtx)
        # delta = delta * (numVtx / 254.0) 

        #basically - move each vertex and update it
        tri_mesh.vertices[uniqueVtx, 0] -= (delta * gradient[uniqueVtx, 0])
        tri_mesh.vertices[uniqueVtx, 1] -= (delta * gradient[uniqueVtx, 1])
        tri_mesh.vertices[uniqueVtx, 2] -= (delta * gradient[uniqueVtx, 2])

        return tri_mesh

    # def meshHFOpt(self, hfObjectiveFcn, calcHFAllMesh, calcMaxHF, calcEnergy, meshObj, coefficientsList, threshold, delta, id):
    def meshHFOpt(self, hfObjectiveFcn, constraint, updateHFProfile, calcHFAllMesh, calcMaxHF, calcEnergy, meshObj, coefficientsList, threshold, delta, fwdModel, id):
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
        vertices = trimeshSolid.vertices
        faces = trimeshSolid.faces
        face_adjacency = trimeshSolid.face_adjacency
        face_adjacency_edges = trimeshSolid.face_adjacency_edges
        initialVolume = trimeshSolid.volume
        print(f"Initial volume: {initialVolume}")
        initialParams = [initialVolume] #[initialVolume, vertices, faces, face_adjacency, face_adjacency_edges]
        # objFcn = hfObjectiveFcn(trimeshSolid, coefficientsList, facesToMove)
        objFcn = hfObjectiveFcn(vertices, faces, face_adjacency, face_adjacency_edges, initialParams, coefficientsList, facesToMove)
        # objFcn = hfObjectiveFcn(trimeshSolid.vertices, trimeshSolid.faces, coefficientsList, facesToMove)
        # all_objective_function_values = [hfObjectiveFcn(trimeshSolid, coefficientsList, facesToMove)]
        all_objective_function_values = [objFcn[0]]
        all_sum_normals_diff = [objFcn[1]]
        all_max_normals_diff = [objFcn[2]]
        hf_all_mesh = calcHFAllMesh(fwdModel.hfMode, fwdModel.q_dir, fwdModel.q_mag, trimeshSolid) #calcHFAllMesh(trimeshSolid)
        max_hf_each_run = [calcMaxHF(hf_all_mesh, facesToMove)]

        # make VTK to display HF on surface
        #self.plotHFVTK(calcHFAllMesh(trimeshSolid), trimeshSolid, f"{id}", count=4242)
        self.plotHFVTK(calcHFAllMesh(fwdModel.hfMode, fwdModel.q_dir, fwdModel.q_mag, trimeshSolid), trimeshSolid, f"{id}", count=4242)

        prev_objVal = 2000
        curr_objVal = 0

        # t0 = time.time()

        #faces to NOT move
        facesToKeep = indicesToNotMove

        while abs(prev_objVal - curr_objVal) > threshold and count < 2000:

            #hf_all_mesh = calcHFAllMesh(trimeshSolid)
            hf_all_mesh = calcHFAllMesh(fwdModel.hfMode, fwdModel.q_dir, fwdModel.q_mag, trimeshSolid)

            t1 = time.time()

            #calc the gradient
            trimeshSolid = self.gradientDescentHF(trimeshSolid, hfObjectiveFcn, hf_all_mesh, face_adjacency, face_adjacency_edges, initialParams, facesToKeep, facesToMove, coefficientsList, delta, f"test{id}", count)

            # print(f"{count}: gradient descent time: {time.time() - t0}")
            print(f"{count}: gradient descent time: {time.time() - t1}")

            vertices = trimeshSolid.vertices

            # input()
            
            #recalculate the hf profile on the surface - don't need this for spheretests
            # updateHFProfile(trimeshSolid) 

            # new_objVal = hfObjectiveFcn(trimeshSolid, coefficientsList, facesToMove)
            # t2 = time.time()
            # new_objVal_all = hfObjectiveFcn(trimeshSolid, coefficientsList, facesToMove)
            new_objVal_all = hfObjectiveFcn(vertices, faces, face_adjacency, face_adjacency_edges, initialParams, coefficientsList, facesToMove)
            # new_objVal_all = hfObjectiveFcn(trimeshSolid.vertices, trimeshSolid.faces, coefficientsList, facesToMove)
            new_objVal = new_objVal_all[0]
            new_sum_normals_diff = new_objVal_all[1]
            new_max_normals_diff = new_objVal_all[2]
            all_objective_function_values.append(new_objVal)
            all_sum_normals_diff.append(new_sum_normals_diff)
            all_max_normals_diff.append(new_max_normals_diff)

            prev_objVal = curr_objVal
            curr_objVal = new_objVal

            # print(f"{count}: time to update all calcs/plots/everything else: {time.time() - t2}")

            # input()

            # new_max_hf = calcMaxHF(hf_all_mesh, facesToMove)
            # max_hf_each_run.append(new_max_hf)

            ##this is for plotting surface fit onto mesh
            # if count % 5 == 0:
            #     self.makePolyFitSurface(trimeshSolid.vertices[constrainedVIdx], f"{id}", count)

            # print(f"New objective function value: {new_objVal}")

            if count and count % 10 == 0: #count % 5 == 0: 
                #self.plotHFVTK(calcHFAllMesh(trimeshSolid), trimeshSolid, f"{id}", count)
                self.plotHFVTK(calcHFAllMesh(fwdModel.hfMode, fwdModel.q_dir, fwdModel.q_mag, trimeshSolid), trimeshSolid, f"{id}", count)

            count += 1
        
        # self.plotRun(all_objective_function_values, max_hf_each_run, f"{id}")
        # self.plotHFVTK(calcHFAllMesh(trimeshSolid), trimeshSolid, f"{id}", count)
        self.plotHFVTK(calcHFAllMesh(fwdModel.hfMode, fwdModel.q_dir, fwdModel.q_mag, trimeshSolid), trimeshSolid, f"{id}", count)
        self.plotMaxNormalsDiff(all_max_normals_diff, f"{id}")
        self.plotNormalsDiff(all_sum_normals_diff, f"{id}")
        self.plotObjectiveFunction(all_objective_function_values, f"{id}")
        # self.plotMaxHF(max_hf_each_run, f"{id}")  
          
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