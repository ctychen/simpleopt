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
        x_allVals = []
        y_allVals = []
        z_allVals = []
        #for each vertex
        for i in range(len(tri_mesh.vertices)):
            #for each dimension
            x_allVals.append(tri_mesh.vertices[i][0])
            y_allVals.append(tri_mesh.vertices[i][1])
            z_allVals.append(tri_mesh.vertices[i][2])
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
        
        # print(f"x values: {x_allVals}")
        # print(f"y values: {y_allVals}")
        # print(f"length x: {len(x_allVals)} length y: {len(y_allVals)} length z: {len(z_allVals)} length gradient: {len(gradient)}")
        return [gradient, x_allVals, y_allVals, z_allVals]
    
    
    def plotIteration(self, gradDescOut, count, fileID):
        # Create a 3D surface plot with color mapped to function values
        gradient = gradDescOut[0]
        x_vals = gradDescOut[1]
        y_vals = gradDescOut[2]
        z_vals = gradDescOut[3]

        gradient_magnitudes = []
        for grad in gradient: 
            gradient_magnitudes.append(np.linalg.norm(grad))

        fig = go.Figure(data=go.Volume(
            x=x_vals,
            y=y_vals,
            z=z_vals,
            value=gradient_magnitudes,
            isomin=np.min(gradient_magnitudes),
            isomax=np.max(gradient_magnitudes),
            opacity=0.3,
            surface_count=17,
            colorscale='Thermal'
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

    def meshHFOpt(self, hfObjectiveFcn, meshObj, changeMeshFcn, threshold, delta, id):
    # def meshHFOpt(self, hfFunction, hfObjectiveFcn, meshObj, threshold, step, id):
        """
        runs optimization process until objective fcn value reaches stopping condition @ minimum
        modifies the mesh based on gradient by applying changeMeshFcn accordingly

        can change changeMeshFcn, hfObjectiveFcn to use different functions
        if we want a different manipulation, or add more stuff to the functions
        """
        #assuming input is already a trimesh, ie. processing the solid was done already
        trimeshSolid = meshObj
        count = 0
        all_objective_function_values = [hfObjectiveFcn(trimeshSolid)]

        print("Starting the mesh HF opt")

        print(f"Starting max HF: {hfObjectiveFcn(trimeshSolid)}")

        #objective fcn should take in a trimesh mesh object
        while hfObjectiveFcn(trimeshSolid) > threshold:

            print(f"Current max HF: {hfObjectiveFcn(trimeshSolid)}")

            #calc the gradient
            gradientDescentOut = self.gradientDescentHF(trimeshSolid, hfObjectiveFcn, delta, fileID=id)
            gradient = gradientDescentOut[0]
            print(f"Gradient calculated: {gradient}")

            self.plotIteration(gradientDescentOut, count, id)

            #move the vertices a bit based on the gradient
            trimeshSolid.vertices = changeMeshFcn(trimeshSolid, gradient, delta)

            trimeshSolid.export(f"test{id}/{count}.stl")

            new_obj_val = hfObjectiveFcn(trimeshSolid)
            all_objective_function_values.append(new_obj_val)

            print(f"New HF value: {new_obj_val}")
            count += 1
        
        #when process is done, the mesh should have been modified - so return it 
        return trimeshSolid
    

    def plotRun(self, objective_function_values, outputDir):
        """
        plot values of objective function over iterations
        """
        x_count = np.linspace(0, len(objective_function_values), len(objective_function_values))
        fig = px.scatter(x = x_count, y = objective_function_values)
        fig.update_xaxes(title_text='Iterations')
        fig.update_yaxes(title_text='Sum of HF over all elements')
        fig.show()            
        output_file = f"{outputDir}/entire_run.html"
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

