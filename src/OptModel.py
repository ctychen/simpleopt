import FreeCAD
import Part
import Mesh
import MeshPart
from FreeCAD import Base
import numpy as np

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import plotly.express as px


class OptModel_MeshHF: 
    """
    Mesh element manipulation to optimize heat flux - model for optimization
    """

    def __init__(self):
        """
        TODO: figure out what properties we actually need for optmodel, and what can just be fcn inputs
        """
        return
    

    def gradientDescentHF(self, tri_mesh, objectiveFunction, delta, fileID):
        """
        gradient descent implementation for heat flux minimization
        takes in trimesh object, applies objective function to it
        moves all mesh vertices by a small amount and finds the gradient when doing so 
        """
        gradient = np.zeros_like(tri_mesh.vertices)
        #for each vertex
        for i in range(len(tri_mesh.vertices)):
            #for each dimension
            for j in range(3):
                #move vertex a bit, see if the objective function decreases
                # print(f"Vertex: {tri_mesh.vertices[i]}") #format: [x, y, z]
                tri_mesh.vertices[i, j] += delta
                # print(f"After moving: {tri_mesh.vertices[i, j]}")
                obj_afterMoving = objectiveFunction(tri_mesh)
                #tri_mesh.export(f"test{fileID}/wip_movingvertices_{i}_{j}.stl")
                #move vertex back to original, calc objective function then for comparison
                tri_mesh.vertices[i, j] -= delta
                obj_beforeMoving = objectiveFunction(tri_mesh)
                gradient[i, j] = (obj_afterMoving - obj_beforeMoving) / (2 * delta)
                # print(f"Gradient value: {gradient[i, j]}")
                # print(f"Gradient should be: {(obj_afterMoving - obj_beforeMoving) / (2 * delta)}")
        
        return gradient

    
    def gradientDesignHFV2(self, tri_mesh, objectiveFunction, allmeshelementsHF, delta, filedir, count):
        """
        what this should do: 
        find all vertices for mesh elements that have heat flux that's maxHF-delta to maxHF
        move those first - calc gradient and move actually depending on gradient
        then do gradient descent on rest

        """ 

        #find mesh elements with HF close to maxHF and calculate those motions first, then apply them
        #modify those within this loop before calculating changes on everything else
        #then set the gradient values for those indices to be 0 afterwards so they don't move while other elements on pass 2 will
    
        #which means calculate hf for all faces using hfCalc (or get it from facesHF / q_all as an input)
        gradient_step0 = np.zeros_like(tri_mesh.vertices)

        highestHFVertexIndices = []

        #maxHF = np.max(allmeshelementsHF)
        
        #changing mesh to refer to element with the most negative HF - maybe this is why the back face is moving?
        maxHF = np.min(allmeshelementsHF)

        print(f"Starting gradient descent. maxHF is {maxHF}")
        
        for i in range(len(allmeshelementsHF)):
            # if allmeshelementsHF[i] > (maxHF - 0.5): #arbitrary threshold we can change later
            if allmeshelementsHF[i] < (maxHF + 0.5):
                i_vertices = tri_mesh.faces[i]  #format: [x, y, z]
                #print(f"from this face: {i_vertices}")

                for v_idx in i_vertices:
                    if v_idx not in highestHFVertexIndices:
                        highestHFVertexIndices.append(v_idx)
                #if value in my_array[:, col_num]

        print(f"indices found with close-to-maxHF values: {highestHFVertexIndices}")

        for idx in highestHFVertexIndices: 
            for j in range(3):
                tri_mesh.vertices[idx, j] += delta
                obj_afterMoving = objectiveFunction(tri_mesh)
                tri_mesh.vertices[idx, j] -= delta
                obj_beforeMoving = objectiveFunction(tri_mesh)
                gradient_step0[idx, j] = (obj_afterMoving - obj_beforeMoving) / (2 * delta)

        print(f"Found gradients for highest HF elements")

        #now, apply gradient to those vertices first so that it is already moved 
        for idx in highestHFVertexIndices: 
            for j in range(3):
                tri_mesh.vertices[idx, j] -= (delta * gradient_step0[idx, j])

        print(f"Moved highest HF elements")

        #make stl after moving highest-hf-elements to check what's happening
        tri_mesh.export(f"{filedir}/{count}_onlyhighesthfmoved.stl")

        gradient = np.zeros_like(tri_mesh.vertices)
        #for each vertex
        for i in range(len(tri_mesh.vertices)):
            #for each dimension
            for j in range(3):
                #now, calc gradient if we haven't already moved the points
                if i not in highestHFVertexIndices: 
                    #move vertex a bit, see if the objective function decreases
                    tri_mesh.vertices[i, j] += delta
                    obj_afterMoving = objectiveFunction(tri_mesh)
                    #move vertex back to original, calc objective function then for comparison
                    tri_mesh.vertices[i, j] -= delta
                    obj_beforeMoving = objectiveFunction(tri_mesh)
                    gradient[i, j] = (obj_afterMoving - obj_beforeMoving) / (2 * delta)
                else: #if we've already moved the vertices, keep this 0 since we don't want to move it more? maybe not... 
                    tri_mesh.vertices[i, j] = 0.0 

            #testing out stuff below this is new
            tri_mesh.vertices[i, 0] -= (delta * gradient[i, 0])
            tri_mesh.vertices[i, 1] -= (delta * gradient[i, 1])
            tri_mesh.vertices[i, 2] -= (delta * gradient[i, 2])

            tri_mesh.export(f"{filedir}/{count}_vertex_{i}_moved.stl")

            input(f"moved vertex {i}")

        #now, move all the vertices that haven't already been moved based on gradient
        return tri_mesh

        #return gradient 


    

    #this doesn't work yet
    def plotIteration(self, gradDescOut, count, fileID):
        # Create a 3D surface plot with color mapped to function values
        gradient = gradDescOut[0]
        x_vals = gradDescOut[1]
        y_vals = gradDescOut[2]
        z_vals = gradDescOut[3]

        gradient_magnitudes = []
        for grad in gradient: 
            print(f"Normal: {np.linalg.norm(grad)}")
            gradient_magnitudes.append(np.linalg.norm(grad))

        # Assuming you have 1D lists of coordinates and magnitudes
        x_coords = [x for x in range(len(x_vals)) for _ in range(len(y_vals)) for __ in range(len(z_vals))]
        y_coords = [y for _ in range(len(x_vals)) for y in range(len(y_vals)) for __ in range(len(z_vals))]
        z_coords = [z for _ in range(len(x_vals)) for __ in range(len(y_vals)) for z in range(len(z_vals))]
        magnitudes = [gradient_magnitudes for _ in range(len(x_vals)) for __ in range(len(y_vals)) for ___ in range(len(z_vals))]

        # Convert the lists to numpy arrays
        x_coords_np = np.array(x_coords)
        y_coords_np = np.array(y_coords)
        z_coords_np = np.array(z_coords)
        magnitudes_np = np.array(magnitudes)

        # Reshape the 1D numpy arrays to 3D
        x_coords_3d = x_coords_np.reshape((len(x_vals), len(y_vals), len(z_vals)))
        y_coords_3d = y_coords_np.reshape((len(x_vals), len(y_vals), len(z_vals)))
        z_coords_3d = z_coords_np.reshape((len(x_vals), len(y_vals), len(z_vals)))
        magnitudes_3d = magnitudes_np.reshape((len(x_vals), len(y_vals), len(z_vals)))


        fig = go.Figure(data=go.Isosurface(
            x=x_coords_3d.flatten(),
            y=y_coords_3d.flatten(),
            z=z_coords_3d.flatten(),
            value=magnitudes_3d.flatten(),
            opacity=0.6,
            isomin=magnitudes_3d.min(),
            isomax=magnitudes_3d.max(),
            surface_count=3,
            caps=dict(x_show=False, y_show=False)
            ))

        # Set plot layout and axis labels
        fig.update_layout(
            title=f'Vertices colored by gradient magnitude value on iteration {count}',
            scene=dict(
                xaxis_title='X',
                yaxis_title='Y',
                zaxis_title='Z'
            )
        )

        # Show the plot
        fig.show()

        output_file = f'test{fileID}/{count}_plot_of_iteration.html'
        pio.write_html(fig, output_file)

        print(f"Plotted iteration")
        return
    

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
        #TODO: 

        #assuming input is already a trimesh, ie. processing the solid was done already
        trimeshSolid = meshObj
        count = 0
        all_objective_function_values = [hfObjectiveFcn(trimeshSolid)]
        max_hf_each_run = [calcMaxHF(trimeshSolid)]
        sum_hf_each_run = [calcHFSum(trimeshSolid)] #[hfAllMesh(trimeshSolid)]

        print("Starting the mesh HF opt")

        print(f"Starting objective function value: {hfObjectiveFcn(trimeshSolid)}")
        trimeshSolid.export(f"test{id}/original.stl")

        prev_objVal = 2000
        curr_objVal = 0

        while abs(prev_objVal - curr_objVal) > threshold: #and prev_objVal > curr_objVal: #or not(constraints(trimeshSolid)): 

        #objective fcn should take in a trimesh mesh object --> below was original threshold
        #while hfObjectiveFcn(trimeshSolid) > threshold:

            #v2: before calc gradient, calc the max hf so that elements with it can be moved first
            #maxHF = calcMaxHF(trimeshSolid)
            hf_all_mesh = calcHFAllMesh(trimeshSolid)
            
            #calc the gradient
            #gradient = self.gradientDesignHFV2(trimeshSolid, hfObjectiveFcn, hf_all_mesh, delta, f"test{id}", count)
            newTrimesh = self.gradientDesignHFV2(trimeshSolid, hfObjectiveFcn, hf_all_mesh, delta, f"test{id}", count)
            # gradient = self.gradientDescentHF(trimeshSolid, hfObjectiveFcn, hf_all_mesh, delta, fileID=id)

            trimeshSolid = newTrimesh

            #move the vertices a bit based on the gradient
            #trimeshSolid.vertices = changeMeshFcn(trimeshSolid, gradient, delta)

            trimeshSolid.export(f"test{id}/{count}.stl")

            new_objVal = hfObjectiveFcn(trimeshSolid)
            all_objective_function_values.append(new_objVal)

            prev_objVal = curr_objVal
            curr_objVal = new_objVal

            new_max_hf = calcMaxHF(trimeshSolid)
            max_hf_each_run.append(new_max_hf)

            new_sum_hf = calcHFSum(trimeshSolid) #hfAllMesh(trimeshSolid)
            sum_hf_each_run.append(new_sum_hf)

            print(f"New objective function value: {new_objVal}")

            # #Plotting process, for anything using means/distributions
            # q_mesh_all = hfAllMesh(trimeshSolid)
            # fig = go.Figure(data=[go.Histogram(x=q_mesh_all)])
            # fig.show()
            # output_file = f"test{id}/{count}_run_entiredistribution.html"
            # pio.write_html(fig, output_file)

            # q_mesh_objective = calcMaxHF(trimeshSolid)
            # fig = go.Figure(data=[go.Histogram(x=q_mesh_objective)])
            # fig.show()
            # output_file = f"test{id}/{count}_run_distributionforobjective.html"
            # pio.write_html(fig, output_file)

            if count and count % 1 == 0: 
                x_count = np.linspace(0, len(all_objective_function_values), len(all_objective_function_values))
                fig = px.scatter(x = x_count, y = all_objective_function_values)
                fig.update_xaxes(title_text='Iterations')
                fig.update_yaxes(title_text='Objective function: c1*max(HF) + c2*(sum HF over elements)')
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


            count += 1
        
        self.plotRun(all_objective_function_values, max_hf_each_run, sum_hf_each_run, f"test{id}")
        
        #when process is done, the mesh should have been modified - so return it 
        return trimeshSolid
    

    def plotRun(self, objective_function_values, max_hf_each_run, sum_hf_each_run, outputDir):
        """
        plot values of objective function over iterations
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

