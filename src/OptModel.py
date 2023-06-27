import FreeCAD
import Part
import Mesh
import MeshPart
from FreeCAD import Base
import numpy as np

import time

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import plotly.express as px


import vtk
from vtk import vtkPolyData, vtkPoints, vtkCellArray, vtkDoubleArray, vtkPolyDataWriter, vtkTriangle

class OptModel_MeshHF: 
    """
    Mesh element manipulation to optimize heat flux - model for optimization
    """

    def __init__(self):
        """
        TODO: figure out what properties we actually need for optmodel, and what can just be fcn inputs
        """
        return
    
    #simulated annealing attempt???

    def simulated_annealing(self, tri_mesh, objective_function, calcHFAllMesh, calcMaxHFMesh, delta, initial_temperature, cooling_rate, threshold, id):
        """
        simulated annealing attempt for heat flux minimization
        """

        # Initialization
        current_solution = tri_mesh
        current_objective = objective_function(current_solution)
        
        temperature = initial_temperature
        count = 0
        all_objectives = [current_objective]
        all_maxHF = [calcMaxHFMesh(current_solution)]

        # Iteration
        while temperature > threshold:  # Stop when the temperature is low
            print(f"Temperature: {temperature}")
            # Generate a neighbor
            neighbor = self.generate_neighbor(current_solution, delta, calcHFAllMesh)

            # Evaluate the neighbor
            neighbor_objective = objective_function(neighbor)

            # self.plotHFVTK(calcHFAllMesh(neighbor), neighbor, f"test{id}/neighbors/", count)

            # Decide whether to move
            if (neighbor_objective < current_objective or 
                np.random.rand() < np.exp(-(neighbor_objective - current_objective) / temperature)):
                current_solution = neighbor
                current_objective = neighbor_objective

                current_maxHF = calcMaxHFMesh(current_solution)

                all_objectives.append(current_objective)
                all_maxHF.append(current_maxHF)

                self.plotHFVTK(calcHFAllMesh(current_solution), current_solution, f"test{id}", count)

            # Update the temperature according to the cooling schedule
            temperature *= cooling_rate

            if count and count % 5 == 0: 
                x_count = np.linspace(0, len(all_objectives), len(all_objectives))
                fig = px.scatter(x = x_count, y = all_objectives)
                fig.update_xaxes(title_text='Iterations')
                fig.update_yaxes(title_text=f'Objective function: {objective_function.__name__}')
                fig.show()            
                output_file = f"test{id}/objective_up_to_run_{count}.html"
                pio.write_html(fig, output_file)

                x_count = np.linspace(0, len(all_maxHF), len(all_maxHF))
                fig = px.scatter(x = x_count, y = all_maxHF)
                fig.update_xaxes(title_text='Iterations')
                fig.update_yaxes(title_text='Max HF')
                fig.show()            
                output_file = f"test{id}/max_hf_up_to_run_{count}.html"
                pio.write_html(fig, output_file)

            count += 1

        return current_solution

    def generate_neighbor(self, tri_mesh, delta, calcAllMeshHF):
        """
        generate neighbor solution by moving the vertices slightly
        """
        neighbor = tri_mesh.copy()

        allmeshelementsHF = calcAllMeshHF(neighbor)

        use_set = set(np.where(allmeshelementsHF > 0.0)[0])
        # Sort indices based on allmeshelementsHF values in descending order
        sortedFaceIndices = np.argsort(allmeshelementsHF)[::-1]

        for idx in sortedFaceIndices: 
            if idx in use_set: 
                face = neighbor.faces[idx]

                for vertexIdx in face:  #vertexIdx is VERTEX INDICES
                    displacement = (np.random.rand(3) - 0.5) * 2 * delta 
                    neighbor.vertices[vertexIdx, 0] += displacement[0]
                    neighbor.vertices[vertexIdx, 1] += displacement[1]
                    neighbor.vertices[vertexIdx, 2] += displacement[2]
        
        # Add a small random displacement to each vertex
        # for vertex in neighbor.vertices:
        #     displacement = (np.random.rand(3) - 0.5) * 2 * delta  # Random value between -delta and delta
        #     #vertex += displacement
        #     vertex[0] += displacement[0]
        #     vertex[1] += displacement[1]
        #     vertex[2] += displacement[2]

        return neighbor

    

    def gradientDescentHF(self, tri_mesh, objectiveFunction, allmeshelementsHF, delta, filedir, count):
        """
        gradient descent implementation for heat flux minimization
        takes in trimesh object and sorts elements by HF to deal with worst elements first
        calc gradient for each element by moving vertices a small amount and finding change in objective function
        move each vertex based on gradient * delta when all gradients calculated
        """ 

        #process now:
        #find all vertices and all heat fluxes for all faces
        #make array of 0s with same length as list of faces
        #sort heat flux list by magnitude
        #assign face (list of vertices) to same index that its hf is at so those stay linked
        #convert list of faces to list of vertices, keeping the order
        #this way the faces with the highest hf's get moved/gradients calculated first 
        #

        use_set = set(np.where(allmeshelementsHF > 0.0)[0])

        # Sort indices based on allmeshelementsHF values in descending order
        sortedFaceIndices = np.argsort(allmeshelementsHF)[::-1]

        gradient = np.zeros_like(tri_mesh.vertices)

        #variable delta
        max_delta = delta
        min_delta = 1e-6
        delta_decrease_factor = 0.5
        delta_increase_factor = 1.01
        successful_moves = 0

        normals_threshold = 0.01
        moves_threshold = 5

        # Original position and normals
        # original_vertices = tri_mesh.vertices.copy()

        for idx in sortedFaceIndices: 
            if idx in use_set: 
                face = tri_mesh.faces[idx]

                for vertexIdx in face:  #vertexIdx is VERTEX INDICES

                    obj_beforeMoving = objectiveFunction(tri_mesh)
                    original_normals = tri_mesh.face_normals
                    # print(f"Original normals: {original_normals}")

                    for j in range(3): #for every dimension - move the vertex a bit and calculate the change in objectiveFunction

                        tri_mesh.vertices[vertexIdx, j] += delta
                        obj_afterMoving = objectiveFunction(tri_mesh)

                        tri_mesh.vertices[vertexIdx, j] -= delta

                        gradient[vertexIdx, j] = (obj_afterMoving - obj_beforeMoving) / (2 * delta)
                
                    #basically - move each vertex and update it

                    tri_mesh.vertices[vertexIdx, 0] -= (delta * gradient[vertexIdx, 0])
                    tri_mesh.vertices[vertexIdx, 1] -= (delta * gradient[vertexIdx, 1])
                    tri_mesh.vertices[vertexIdx, 2] -= (delta * gradient[vertexIdx, 2])    

                    #everything below is for checking if changes in normals is within amount we're ok with 
                    proposed_normals = tri_mesh.face_normals
                    # print(f"New normals: {proposed_normals}")
                    
                    max_change = np.max(np.abs(proposed_normals - original_normals))
                    # print(f"Max change in normals: {max_change}")

                    if max_change > normals_threshold:
                        #if mmax change in normals is over the threshold, reject move, decrease delta
                        tri_mesh.vertices[vertexIdx, 0] += (delta * gradient[vertexIdx, 0])
                        tri_mesh.vertices[vertexIdx, 1] += (delta * gradient[vertexIdx, 1])
                        tri_mesh.vertices[vertexIdx, 2] += (delta * gradient[vertexIdx, 2]) 
                        delta = max(delta * delta_decrease_factor, min_delta)
                        print(f"Delta in normals was too high: {max_change}, need to decrease delta to: {delta}")
                        successful_moves = 0  #reset counter for good moves 
                    else:
                        #if max change in normals is below the threshold, accept move, potentially increase delta
                        original_normals = proposed_normals
                        successful_moves += 1
                        if successful_moves >= moves_threshold:  #increase delta after several successful moves
                            delta = min(delta * delta_increase_factor, max_delta)
                            print(f"Hit enough successful moves and increasing delta to: {delta}")
                            successful_moves = 0  #reset counter for good moves

        return tri_mesh


    def moveMeshVertices(self, trimeshSolid, gradient, delta):
        """
        function for how we want to adjust mesh vertices, depending on what the gradient is 
        """
        return trimeshSolid.vertices - (delta * gradient)


    def meshHFOpt(self, hfObjectiveFcn, calcHFAllMesh, calcMaxHF, calcHFSum, meshObj, changeMeshFcn, threshold, delta, id):
    # def meshHFOpt(self, hfFunction, hfObjectiveFcn, meshObj, threshold, step, id):
        """
        runs optimization process until objective fcn value reaches stopping condition @ minimum
        modifies the mesh based on gradient by applying changeMeshFcn accordingly

        can change changeMeshFcn, hfObjectiveFcn to use different functions
        if we want a different manipulation, or add more stuff to the functions
        """
        #TODO: add constraints somehow - take in list of criteria? eg. don't move face if x=0 or y=0 or x=10 or y=10?

        #assuming input is already a trimesh, ie. processing the solid was done already
        trimeshSolid = meshObj

        count = 0

        all_objective_function_values = [hfObjectiveFcn(trimeshSolid)]
        max_hf_each_run = [calcMaxHF(trimeshSolid)]
        sum_hf_each_run = [calcHFSum(trimeshSolid)] 

        print("Starting the mesh HF opt")

        print(f"Starting objective function value: {hfObjectiveFcn(trimeshSolid)}")
        #trimeshSolid.export(f"test{id}/original.stl")

        prev_objVal = 2000
        curr_objVal = 0

        t0 = time.time()

        while abs(prev_objVal - curr_objVal) > threshold: 

            hf_all_mesh = calcHFAllMesh(trimeshSolid)
            
            #calc the gradient
            newTrimesh = self.gradientDescentHF(trimeshSolid, hfObjectiveFcn, hf_all_mesh, delta, f"test{id}", count)

            print(f"Time elapsed for GD {count}: {time.time() - t0}")

            trimeshSolid = newTrimesh

            # trimeshSolid.export(f"test{id}/{count}.stl")

            new_objVal = hfObjectiveFcn(trimeshSolid)
            all_objective_function_values.append(new_objVal)

            prev_objVal = curr_objVal
            curr_objVal = new_objVal

            new_max_hf = calcMaxHF(trimeshSolid)
            max_hf_each_run.append(new_max_hf)

            new_sum_hf = calcHFSum(trimeshSolid) #hfAllMesh(trimeshSolid)
            sum_hf_each_run.append(new_sum_hf)

            #make VTK to display HF on surface
            self.plotHFVTK(calcHFAllMesh(trimeshSolid), trimeshSolid, f"test{id}", count)

            print(f"New objective function value: {new_objVal}")

            if count and count % 5 == 0: 
                x_count = np.linspace(0, len(all_objective_function_values), len(all_objective_function_values))
                fig = px.scatter(x = x_count, y = all_objective_function_values)
                fig.update_xaxes(title_text='Iterations')
                fig.update_yaxes(title_text=f'Objective function: {hfObjectiveFcn.__name__}')
                fig.show()            
                output_file = f"test{id}/objective_up_to_run_{count}.html"
                pio.write_html(fig, output_file)

                x_count = np.linspace(0, len(max_hf_each_run), len(max_hf_each_run))
                fig = px.scatter(x = x_count, y = max_hf_each_run)
                fig.update_xaxes(title_text='Iterations')
                fig.update_yaxes(title_text='Max HF')
                fig.show()            
                output_file = f"test{id}/max_hf_up_to_run_{count}.html"
                pio.write_html(fig, output_file)

                x_count = np.linspace(0, len(sum_hf_each_run), len(sum_hf_each_run))
                fig = px.scatter(x = x_count, y = sum_hf_each_run)
                fig.update_xaxes(title_text='Iterations')
                fig.update_yaxes(title_text='Sum HF on mesh')
                fig.show()            
                output_file = f"test{id}/sum_hf_up_to_run_{count}.html"
                pio.write_html(fig, output_file)

                # #make VTK to display HF on surface
                # self.plotHFVTK(calcHFAllMesh(trimeshSolid), trimeshSolid, f"test{id}")

            count += 1
        
        self.plotRun(all_objective_function_values, max_hf_each_run, sum_hf_each_run, f"test{id}")
        
        #when process is done, the mesh should have been modified - so return it 
        return trimeshSolid


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
        writer.SetFileName(f"{fileDir}/{count}_hfOnMesh.vtk")
        writer.SetInputData(polydata)
        writer.Write()

        return 
    

    def plotRun(self, objective_function_values, max_hf_each_run, sum_hf_each_run, outputDir):
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

        x_count = np.linspace(0, len(sum_hf_each_run), len(sum_hf_each_run))
        fig = px.scatter(x = x_count, y = sum_hf_each_run)
        fig.update_xaxes(title_text='Iterations')
        fig.update_yaxes(title_text='Sum HF on mesh')
        fig.show()            
        output_file = f"{outputDir}/sum_hf_each_run.html"
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

