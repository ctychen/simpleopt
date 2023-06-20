import CADClass
import FreeCAD 
import Part
import Mesh
import MeshPart
import os
import numpy as np
import toolsClass
import math 

import trimesh
# Check available backends
print(trimesh.interfaces.blender.exists)
print(trimesh.interfaces.scad.exists)

tools = toolsClass.tools()

class Box_Vector_Mesh(CADClass.CAD):

    def __init__(self, stlfile="", stpfile="", meshres=1.0):
        super(CADClass.CAD, self).__init__()
        self.STLfile = stlfile
        self.STPfile = stpfile
        self.parts = self.loadSTEP()

        self.meshres = meshres

        self.allmeshes = self.part2mesh(self.CADparts, meshres)
        # self.allmeshes = self.part2meshStandard(self.CADparts)

        self.vertices = np.array(self.allmeshes[0].Points)
        self.faces = np.array(self.allmeshes[0].Facets)

        return

    def makeMesh(self, vertices):
        """
        create freecad mesh from array or list of vertices
        """
        if type(vertices) != list:
            vertices = vertices.tolist()

        newMesh = Mesh.Mesh(vertices)
        return newMesh

    def isPointOutsidePyramid(self, point, pyramidVertices):
        """
        check if given vertex is outside the pyramid defined by its vertices
        not needed if using trimesh implementation
        """
        base = pyramidVertices[0:4]
        apex = pyramidVertices[4]
        
        for i in range(4):
            #vertices for pyramid edge
            v1, v2 = base[i], base[(i+1)%4]
            
            #find normal between edges
            normal = np.cross(v2 - v1, apex - v1)
            
            #vector from the point to the triangle face
            vec = v1 - point
            
            #if dot product is positive, the point is outside pyramid
            if np.dot(normal, vec) >= 0:
                return False
        
        return True
        
    def pointTriangleDistance(self, P, A, B, C):
        """
        calculate distance from a point to a triangle's center using barycentric coordinates
        triangle defined by its vertices
        this was used for finding the distance from a mesh vertex to the nearest surface on the target pyramid
        deprecated if using the trimesh implementation
        """
        # Compute vectors        
        v0 = C - A
        v1 = B - A
        v2 = P - A

        # Compute dot products
        dot00 = np.dot(v0, v0)
        dot01 = np.dot(v0, v1)
        dot02 = np.dot(v0, v2)
        dot11 = np.dot(v1, v1)
        dot12 = np.dot(v1, v2)

        # Compute barycentric coordinates
        invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
        u = (dot11 * dot02 - dot01 * dot12) * invDenom
        v = (dot00 * dot12 - dot01 * dot02) * invDenom

        # Check if point is in triangle
        if (u >= 0) and (v >= 0) and (u + v <= 1):
            proj = A + u*v0 + v*v1
            return np.linalg.norm(P - proj), proj
        else:
            # Here, compute the distance to the edge or vertex
            edge_points = [
                (np.dot(P - A, B - A) / np.dot(B - A, B - A)) * (B - A) + A,
                (np.dot(P - B, C - B) / np.dot(C - B, C - B)) * (C - B) + B,
                (np.dot(P - C, A - C) / np.dot(A - C, A - C)) * (A - C) + C
            ]
            edge_dists = [np.linalg.norm(P - pt) for pt in edge_points]
            min_dist_idx = np.argmin(edge_dists)
            return edge_dists[min_dist_idx], edge_points[min_dist_idx]

    def pointDistanceToPyramid(self, P, pyramidVertices):
        """
        calculate the nearest distance from some point to the pyramid, as defined by its vertices. 
        eventually this was replaced by trimesh nearest_on_surface, but used in non-trimesh implementation.
        """
        base = pyramidVertices[:4]
        apex = pyramidVertices[4]
        
        results = [
            self.pointTriangleDistance(P, base[i], base[(i+1)%4], apex)
            for i in range(4)
        ]
        
        results.append(self.pointTriangleDistance(P, base[0], base[1], base[2]))
        results.append(self.pointTriangleDistance(P, base[2], base[3], base[0]))
        
        min_distance, closest_point = min(results, key=lambda x: x[0])
        
        return min_distance, closest_point
    
    def gradMoveVertex(self, vertex, pyramidVertices, step_size=0.1):
        """
        gradient descent implementation, given pyramid vertices - not using trimesh
        gradually move a vertex in the direction such that it approaches the closest point on the target mesh
        """
        distance, closest_point = self.pointDistanceToPyramid(vertex, pyramidVertices)
        #calculate direction vector from the vertex to the closest point on the pyramid
        #then normalize
        direction_vector = closest_point - vertex
        # norms = np.linalg.norm(direction_vector)
        # if norms != 0.0: #so we don't have division by 0 issues
        #     direction_vector /= np.linalg.norm(direction_vector)
        #move vertex towards closest point using direction vector
        vertex += step_size * direction_vector
        return vertex
    
    def gradMoveVertex_TriMesh(self, vertex, tri_mesh, step_size = 0.1):
        """
        gradient descent implementation for trimesh mesh
        gradually move a vertex in the direction such that it approaches the closest point on the target mesh
        """
        # closest_point, distance, triangle_id = mesh.nearest.on_surface([point])
        closest_point, distance, triangle_id = tri_mesh.nearest.on_surface([vertex])
        # print(f"Closest point on the mesh: {closest_point[0]}")
        # print(f"Distance from point to mesh: {distance[0]}")
        direction_vector = closest_point[0] - vertex
        #if norms != 0.0: #so we don't have division by 0 issues
        #    direction_vector /= np.linalg.norm(direction_vector)
        #move vertex towards closest point using direction vector
        vertex += step_size * direction_vector
        return vertex
    
    

    #find boolean between mesh and target pyramid, and get the volume of the boolean
    #use this volume as objective function for optimization: 
    #want the volume of the final mesh to be as close to the target volume as possible (diff = 0)
    def pyramidFromCubeV2(self, id='013'):

        print(f"Starting pyramid attempt")

        vertices = []

        for i in range(len(self.allmeshes[0].Facets)):
            facet_points = self.allmeshes[0].Facets[i].Points
            for point in facet_points:
                vertices.append(list(point))      

        vertices = np.array(vertices) #cube vertices now defined

        print(f"Number of vertices: {len(vertices)}")

        os.makedirs(f"pyramidtest{id}")
        meshUpdated = self.makeMesh(vertices)
        self.saveMeshSTL(meshUpdated, f"pyramidtest{id}/before_pyramid", self.meshres)

        faces = np.array([[facet.PointIndices[i] for i in range(3)] for facet in self.faces])
        print(f"Number of faces: {faces.shape}")   

        #approach: 
        #objective to minimize: (volume of mesh outside of target volume) + (volume of target outside of mesh), with intersect? 
        #determine pyramid volume based on defined vertices, & determine where this sits in space
        #how to calculate mesh volume outside of another volume? 
        #somehow move mesh points to reduce the volume: probably similar to grad desc? (how to determine though - move vertex to nearest point inside the mesh?)
        #eg: calculate the distance from a mesh vertex to the closest point inside the target area, then move it there (or just move to some distance within it)
        
        #for vertex in mesh: determine if the vertex is inside/outside target pyramid volume
        #either: minimize volume outside the target volume AND/OR: minimize number of vertex points outside the target volume (ie, want it to get to 0) - could be more straightforward?
        #keep moving the vertices until the number outside the pyramid target is 0

        #defining target pyramid object - then use this to make a solid, then mesh, for reference

        pVertices = [
            FreeCAD.Vector(0,0,0), #corners
            FreeCAD.Vector(10,0,0),
            FreeCAD.Vector(10,10,0),
            FreeCAD.Vector(0,10,0),
            FreeCAD.Vector(5,5,10),  # the apex
        ]
        # Create the base
        pBase = Part.makePolygon(pVertices[:4] + [pVertices[0]])  # makePolygon requires a closed loop, hence append the first vertex at the end

        # Create the faces
        pFaces = [Part.Face(pBase)]
        for i in range(4):
            triangle = Part.makePolygon([pVertices[i], pVertices[(i+1)%4], pVertices[4], pVertices[i]])
            pFaces.append(Part.Face(triangle))

        # Create a shell from the faces for the pyramid
        shell = Part.makeShell(pFaces)

        # Create a solid from the shell for the pyramid
        pyramidSolid = Part.makeSolid(shell)
        
        psolid = FreeCAD.ActiveDocument.addObject("Part::Feature", "Pyramid")
        psolid.Shape = pyramidSolid

        pyramidMesh = self.part2meshStandard(psolid)

        print(f'Pyramid Mesh vertices: {pyramidMesh[0].Points}')
        print(f'Pyramid Mesh faces: {pyramidMesh[0].Facets}')

        pyramidMeshVertices = []
        pyramidMeshFaces = []

        pyramidMeshVertices = np.array([[v.x, v.y, v.z] for v in pyramidMesh[0].Points])
        pyramidMeshFaces = np.array([[f.PointIndices[0], f.PointIndices[1], f.PointIndices[2]] for f in pyramidMesh[0].Facets])

        cubeVertices = np.array([[v.x, v.y, v.z] for v in self.allmeshes[0].Points])
        cubeFaces = np.array([[f.PointIndices[0], f.PointIndices[1], f.PointIndices[2]] for f in self.allmeshes[0].Facets])

        count = 0

        #create TriMesh mesh object of the original cube
        cubeTriMesh = trimesh.Trimesh(
            vertices=cubeVertices, faces=cubeFaces
        )

        cubeTriMesh.export(f"pyramidtest{id}/basecube.stl")

        #create TriMesh mesh object of the target pyramid; we had to do it originally from FreeCAD
        pyramidTriMesh = trimesh.Trimesh(
            vertices=pyramidMeshVertices, faces=pyramidMeshFaces
        )

        pyramidTriMesh.export(f"pyramidtest{id}/basepyramid.stl")

        print(f"Pyramid mesh faces: {pyramidMeshFaces}")

        #in theory what this should do:
        #calculate the volume that's inside the cube/object but outside the target volume
        #calculate the volume that's inside the target voluem but outside the cube/object
        #sum the volume
        #this sum is what should be minimized
        #so it's time for mesh ops - boolean/intersection - and FC may not be the way to go

        #both of these inputs should be trimesh objects
        #objective should be: (volume of mesh) - (volume of intersection of/boolean mesh & pyramid) -> minimize
        #if mesh fits perfectly inside pyramid, then Vmesh = Vpyr, so Vbool = Vmesh
        #so objective should be abs(Vmesh - Vbool(mesh,pyr)) since we want that diff to be close to 0
        #trimesh also has fcn for trimesh.boolean.difference
        #trimesh.boolean.intersection(meshes, engine=None, **kwargs)
        def objectiveFunction(cubeMesh, pyramidMesh):
            intersectionMesh = trimesh.boolean.intersection([cubeMesh, pyramidMesh], engine='scad')
            vIntersect = intersectionMesh.volume
            vMesh = cubeMesh.volume
            volumeDiff = np.abs(vMesh - vIntersect)
            return volumeDiff
        
        #finds sum of distances of cube mesh vertices to the pyramid mesh. ideally this should be as close to 0 as possible
        #this objective function ended up unused, but could be useful later?
        def objectiveFunction2(cubeMesh, pyramidMesh):
            sumDistances = 0.0
            for cubeVertex in cubeMesh.Vertices: 
                closest_point, distance, triangle_id = pyramidMesh.nearest.on_surface([cubeVertex])
                sumDistances += distance[0]
            print(f"Sum of vertex-to-pyramid distances: {sumDistances}")
            return sumDistances
        
        #original: used only volume obj fcn, test000-013
        while objectiveFunction(cubeTriMesh, pyramidTriMesh) > 0:

        #new attempt: use 2 objective functions: distance, and volume, to get rid of the "spike" points
       #while (objectiveFunction(cubeTriMesh, pyramidTriMesh) > 0) and (objectiveFunction2(cubeTriMesh, pyramidTriMesh) > 0):
            #move mesh vertices somehow
            #update the mesh somehow
            
            for i in range(len(cubeVertices)):
                #keep step size fixed for now but move vertex a bit, according to direction/distance from pyramid
                print(f"Original cube vertex {i}: {cubeVertices[i]}")

                # cubeVertices[i] = self.gradMoveVertex(cubeVertices[i], pyramidVertices)
                cubeVertices[i] = self.gradMoveVertex_TriMesh(cubeVertices[i], pyramidTriMesh)
                cubeTriMesh.vertices[i] = cubeVertices[i] #Not sure if this is needed but have to check

                print(f"Modified cube vertex {i}: {cubeVertices[i]}")
                print(f"Modified cube vertex in Trimesh {i}: {cubeTriMesh.vertices[i]}")

                if i == 1243: #253 for 2mm, 1243 for 1mm
                    cubeTriMesh.export(f"pyramidtest{id}/wip_{count}.stl")
                    count += 1

        #export and save the final mesh after modification
        mesh2 = self.makeMesh(vertices)
        self.saveMeshSTL(mesh2, f"pyramidtest{id}/pyramid_test_final_{id}", self.meshres)

        return 

    def gradientDescentHF(self, tri_mesh, objectiveFunction, delta):
        gradient = np.zeros_like(tri_mesh.Vertices)
        #for each vertex
        for i in range(len(tri_mesh.Vertices)):
            #for each dimension
            for j in range(3):
                #move vertex a bit, see if the objective function decreases
                tri_mesh.Vertices[i, j] += delta
                obj_afterMoving = objectiveFunction(tri_mesh)
                #move vertex back to original, calc objective function then for comparison
                tri_mesh.Vertices[i, j] -= delta
                obj_beforeMoving = objectiveFunction(tri_mesh)
                gradient[i, j] = (obj_afterMoving - obj_beforeMoving) / (2 * delta)
        return gradient

    def meshHFOpt(self, id='000'):

        print(f"Starting pyramid attempt")

        vertices = []

        for i in range(len(self.allmeshes[0].Facets)):
            facet_points = self.allmeshes[0].Facets[i].Points
            for point in facet_points:
                vertices.append(list(point))      

        vertices = np.array(vertices) #cube vertices now defined

        print(f"Number of vertices: {len(vertices)}")

        os.makedirs(f"pyramidtest{id}")
        meshUpdated = self.makeMesh(vertices)
        self.saveMeshSTL(meshUpdated, f"pyramidtest{id}/before_pyramid", self.meshres)

        faces = np.array([[facet.PointIndices[i] for i in range(3)] for facet in self.faces])
        # print(f"Original faces: {faces}")    
        print(f"Number of faces: {faces.shape}")   

        cubeVertices = np.array([[v.x, v.y, v.z] for v in self.allmeshes[0].Points])
        cubeFaces = np.array([[f.PointIndices[0], f.PointIndices[1], f.PointIndices[2]] for f in self.allmeshes[0].Facets])
        count = 0

        cubeTriMesh = trimesh.Trimesh(
            vertices=cubeVertices, faces=cubeFaces
        )

        cubeTriMesh.export(f"pyramidtest{id}/basecube.stl")
        
        def hfObjectiveFunction(cubeTrimesh):
            #TODO: techically already exists elsewhere, but gotta trimeshify
            #Calculate normal vector for every mesh facet
            #Calculate (dot product of normal vector, direction of HF)*(magnitude of HF)
            #Store all these in a list, then find the max
            meshAllHF = []
            return meshMaxHF
        

        while hfObjectiveFunction(cubeTriMesh) > 0:
            #calculated the objective function value and it means HF is still not @ min
            #move each vertex in direction to minimize HF - how to determine this/pick direction?  

            #calc the gradient
            gradient = self.gradientDescentHF()

            #move the vertices a bit based on the gradient (not sure if you can do this without looping)
            cubeTriMesh.Vertices -= 0.01 * gradient
        
        mesh2 = self.makeMesh(vertices)

        self.saveMeshSTL(mesh2, f"pyramidtest{id}/pyramid_test_final_{id}", self.meshres)

        return 
    

    def saveMeshSTL(self, mesh, label, resolution, path='./'):
        """
        Writes a mesh object to STL file named by part number.
        If mesh is a list of mesh objects, then write a separate file for
        each mesh object in the list.  Clobbers if overWriteMask is True
        """
        #Check if this is a single file or list and make it a list
        if type(mesh) != list:
            mesh = [mesh]
        if type(label)!= np.ndarray:
            if type(label) != list:
                label = [label]
        if type(resolution) != list:
            resolution=[resolution]*len(mesh)

        #Recursively make dirs for STLs
        print("making STL directory")
        tools.makeDir(path, clobberFlag=False, mode=0o774, UID=-1, GID=-1)

        for i in range(len(mesh)):
            # ___ (3 underdashes) is the str we use to separate mesh name from resolution
            # this MATTERS when we read back in a mesh (see self.loadROIMesh and self.loadIntersectMesh)

            #standard meshing algorithm
            stdList = ['standard', 'Standard', 'STANDARD']
            if resolution[i] in stdList:
                filename = path + label[i] + "___"+resolution[i]+".stl"
            #mefisto meshing algorithm
            else:
                filename = path + label[i] + "___{:.2f}mm.stl".format(float(resolution[i]))
            if os.path.exists(filename): # and self.overWriteMask == False:
                print("Not clobbering mesh file...")
            else:
                print("Writing mesh file: " + filename)
                mesh[i].write(filename)
                os.chmod(filename, 0o774)
                # os.chown(filename, self.UID, self.GID)

        print("\nWrote meshes to files")
        return
    
    