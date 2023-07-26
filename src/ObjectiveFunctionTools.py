import numpy as np

import Solid
import ForwardModel
import trimesh

import pandas as pd
import multiprocessing
from multiprocessing import Pool
from concurrent.futures import ThreadPoolExecutor

class ObjectiveFunctionTools: 

    def __init__(self, chmod=0o774):
        #speed benchmarking tools
        self.testN = 0
        self.testTime = 0
        self.chmod = chmod
        return
    
    def setMeshAndGrids(self, trimeshSolid):
        self.all_faces = trimeshSolid.faces
        self.face_adjacency = trimeshSolid.face_adjacency
        self.face_adjacency_edges = trimeshSolid.face_adjacency_edges
        vertices = trimeshSolid.vertices
        self.currentVerticesGrid = np.array([np.tile(vertices[np.newaxis, :], (len(vertices), 1, 1)) for _ in range(3)])
        self.newVerticesGrid = np.array([np.tile(vertices[np.newaxis, :], (len(vertices), 1, 1)) for _ in range(3)])
        self.initialIMC = self.calculateNonNormalizedIntegralMeanCurvature(vertices, self.all_faces, self.face_adjacency, self.face_adjacency_edges)
        return 
    
    def setParams(self, initialParams, coefficientsList):
        self.initialParams = initialParams
        self.coefficientsList = coefficientsList
        return 
    
    def setFacesToMove(self, facesToMove):
        self.facesToMove = facesToMove
        return 
    
    # def setNewVerticesGrid(self, newVerticesGrid, delta):
    #     self.newVerticesGrid = newVerticesGrid
    #     range_indices = np.arange(self.currentVerticesGrid.shape[1])
    #     self.newVerticesGrid[0, range_indices, range_indices, 0] += delta
    #     self.newVerticesGrid[1, range_indices, range_indices, 1] += delta
    #     self.newVerticesGrid[2, range_indices, range_indices, 2] += delta
    #     return 

    def setNewVerticesGrid(self, newVerticesGrid):
        self.newVerticesGrid = newVerticesGrid
        return 
    
    def calculateFaceNormals(self, vertices, faces):
        """
        calculate normal vectors for each face in a trimesh.

        Parameters
        ----------
        vertices : (n, 3) array
            Vertex coordinate in the space.
        faces : (m, 3) array
            Indices of vertices that make up each face.
        
        Returns
        -------
        face_normals : (m, 3) array
            Normal vector for each face.
        """

        # #get vectors of the two edges for each face
        # vec1 = vertices[faces[:, 1]] - vertices[faces[:, 0]]
        # vec2 = vertices[faces[:, 2]] - vertices[faces[:, 0]]
        
        # #calculate normal vectors using cross product
        # faceNormals = np.cross(vec1, vec2)
        
        # # Normalize each normal vector
        # norms = np.linalg.norm(faceNormals, axis=1)
        # faceNormals = faceNormals / norms[:, np.newaxis]

        #get vectors of the two edges for each face
        vecs = vertices[faces[:, 1:]] - vertices[faces[:, 0, np.newaxis]]
        
        #calculate normal vectors using cross product
        faceNormals = np.cross(vecs[:, 0], vecs[:, 1])
        
        # Normalize each normal vector
        norms = np.linalg.norm(faceNormals, axis=1, keepdims=True)
        faceNormals /= norms
        
        return faceNormals
    
    ###OBJECTIVE FUNCTIONS AND CONSTRAINTS###

    # def objective_for_vertex_dim(objectiveFunction, newVerticesGrid, vtx, dim, all_faces, face_adjacency, face_adjacency_edges, initialParams, coefficientsList, facesToMove):
    #     #face_adjacency, face_adjacency_edges, initialParams
    #     return objectiveFunction(newVerticesGrid[dim][vtx], all_faces, face_adjacency, face_adjacency_edges, initialParams, coefficientsList, facesToMove)[0]
    
    # def objective_function_for_vertex_dim(self, vtx, dim):
        #return self.objectiveFunction(self.newVerticesGrid[dim][vtx], self.all_faces, self.face_adjacency, self.face_adjacency_edges, self.initialParams, self.coefficientsList, self.facesToMove)[0]

    def calculateNormalsDiff(self, trimeshSolid):
        vertex_defects = trimesh.curvature.vertex_defects(trimeshSolid)
        sumVertexDefects = np.sum(np.abs(vertex_defects))
        maxVertexDefect = np.max(np.abs(vertex_defects))  
        maxAngleBetweenNormals = 0 #TODO - will need to bring in these terms from impelmentation in RunModel
        return sumVertexDefects, maxVertexDefect, maxAngleBetweenNormals
    
    def calculateAngles(self, pairs):
        pairs = np.asanyarray(pairs, dtype=np.float64)
        # do the dot product between vectors
        dots = np.dot(pairs[:, 0] * pairs[:, 1], [1.0] * pairs[:, 0].shape[1])
        # clip for floating point error
        dots = np.clip(dots, -1.0, 1.0)
        # do cos and remove arbitrary sign
        angles = np.abs(np.arccos(dots))
        return angles

    def calculateDotProducts(self, vecA, vecB):
        #from trimesh - util.py
        vecA = np.asanyarray(vecA)
        return np.dot(vecA * vecB, [1.0] * vecA.shape[1])
        
    def calculateNonNormalizedIntegralMeanCurvature(self, vertices, faces, face_adjacency, face_adjacency_edges):
        """
        Compute the integral mean curvature of a mesh. Not normalized. This was the original function used
        """
        #calc face normals
        #face_normals = Solid.calculateFaceNormals(vertices, faces)
        face_normals = self.calculateFaceNormals(vertices, faces)
        #calc angles between adjacent faces
        pairs = face_normals[face_adjacency]
        angles = self.calculateAngles(pairs)
        #calc integral mean curvature
        edges_length = np.linalg.norm(np.subtract(
            *vertices[face_adjacency_edges.T]), axis=1)
        integralMeanCurvature = (angles * edges_length).sum() * 0.5
        return integralMeanCurvature
    
    def calculateIntegralMeanCurvature(self, vertices, faces, face_adjacency, face_adjacency_edges, initialIMC):
        """
        Compute the integral mean curvature of a mesh. Normalized: sum(curvatures / sum of curvatures)
        """
        #calc face normals
        #face_normals = Solid.calculateFaceNormals(vertices, faces)
        face_normals = self.calculateFaceNormals(self, vertices, faces)
        #calc angles between adjacent faces
        pairs = face_normals[face_adjacency]
        angles = self.calculateAngles(pairs)
        #calc integral mean curvature
        edges_length = np.linalg.norm(np.subtract(
            *vertices[face_adjacency_edges.T]), axis=1)
        
        allValues = (angles * edges_length) * 0.5
        imc = np.sum(allValues)
        integralMeanCurvature = imc / initialIMC
        #integralMeanCurvature_all = (angles * edges_length).sum() * 0.5
        return integralMeanCurvature   
    
    def calculateIntegralMeanCurvatureParallel(self, vtx, dim):
        """
        Compute the integral mean curvature of a mesh. Normalized: sum(curvatures / sum of curvatures)
        """
        #calc face normals
        vertices = self.newVerticesGrid[dim][vtx]
        #face_normals = Solid.calculateFaceNormals(vertices, self.all_faces)
        face_normals = self.calculateFaceNormals(vertices, self.all_faces)
        #calc angles between adjacent faces
        pairs = face_normals[self.face_adjacency]
        angles = self.calculateAngles(pairs)
        #calc integral mean curvature
        edges_length = np.linalg.norm(np.subtract(
            *vertices[self.face_adjacency_edges.T]), axis=1)
        
        allValues = (angles * edges_length) * 0.5
        imc = np.sum(allValues)
        integralMeanCurvature = imc / self.initialIMC
        return integralMeanCurvature   

    def calcUnitVectors(self, vectors):
        norm = np.sqrt(np.dot(vectors, vectors))
        unitVectors = vectors / norm
        return unitVectors

    def calculateVertexDefects(self, vertices, faces, face_adjacency):
        #calculate mesh face angles - so angles of triangles
        triangles = np.asanyarray(vertices.view(np.ndarray)[faces], dtype=np.float64)
        # get a unit vector for each edge of the triangle
        u = self.calcUnitVectors(triangles[:, 1] - triangles[:, 0])
        v = self.calcUnitVectors(triangles[:, 2] - triangles[:, 0])
        w = self.calcUnitVectors(triangles[:, 2] - triangles[:, 1])

        # run the cosine and per-row dot product
        angles = np.zeros((len(triangles), 3), dtype=np.float64)
        # clip to make sure we don't float error past 1.0
        angles[:, 0] = np.arccos(np.clip(self.calculateDotProducts(u, v), -1, 1))
        angles[:, 1] = np.arccos(np.clip(self.calculateDotProducts(-u, w), -1, 1))
        # the third angle is just the remaining
        angles[:, 2] = np.pi - angles[:, 0] - angles[:, 1]

        # a triangle with any zero angles is degenerate
        # so set all of the angles to zero in that case
        angles[(angles < 1e-8).any(axis=1), :] = 0.0

        angle_sum = np.array(angles.sum(axis=1)).flatten()

        return (2*np.pi) - angle_sum
    

    def vtxFacesObjectiveFunctionCalc(self, verticesList):
        c0, c1, c2, c3, c4 = self.coefficientsList
        imcTerm = c2 * self.calculateIntegralMeanCurvature(verticesList, self.all_faces, self.face_adjacency, self.face_adjacency_edges, self.initialIMC)
        vertexDefectsTerm = 0#c2 * calculateVertexDefects(vertices, faces, face_adjacency)
        maxNormalsTerm = 0#c3 * maxAngleBetweenNormals   
        return imcTerm + maxNormalsTerm + vertexDefectsTerm #[imcTerm + maxNormalsTerm + vertexDefectsTerm, imcTerm, maxNormalsTerm]

    def objectiveFunction(self, vtx, dim): 
        c0, c1, c2, c3, c4 = self.coefficientsList
        # imcTerm = c2 * self.calculateIntegralMeanCurvature(self.newVerticesGrid[dim][vtx], self.all_faces, self.face_adjacency, self.face_adjacency_edges, self.initialIMC)
        imcTerm = c2 * self.calculateIntegralMeanCurvatureParallel(vtx, dim)
        vertexDefectsTerm = 0#c2 * calculateVertexDefects(vertices, faces, face_adjacency)
        maxNormalsTerm = 0#c3 * maxAngleBetweenNormals   
        return imcTerm + maxNormalsTerm + vertexDefectsTerm #[imcTerm + maxNormalsTerm + vertexDefectsTerm, imcTerm, maxNormalsTerm]
    

    # def objective_for_vertex_dim(objectiveFunction, newVerticesGrid, vtx, dim, all_faces, face_adjacency, face_adjacency_edges, initialParams, coefficientsList, facesToMove):
    #     #face_adjacency, face_adjacency_edges, initialParams
    #     return objectiveFunction(newVerticesGrid[dim][vtx], all_faces, face_adjacency, face_adjacency_edges, initialParams, coefficientsList, facesToMove)[0]
    
    # def objective_function_for_vertex_dim(self, vtx, dim):
        #return self.objectiveFunction(self.newVerticesGrid[dim][vtx], self.all_faces, self.face_adjacency, self.face_adjacency_edges, self.initialParams, self.coefficientsList, self.facesToMove)[0]

    #face_adjacency, faces, and face_adjacency_edges are from trimesh but all will only need to be accessed once at beginning, and are all np arrays
    # def objectiveFunction(vertices, faces, face_adjacency, face_adjacency_edges, initialParams, coefficientsList, unconstrainedFaces):
    #     #for initialParams - this is volume of original mesh, etc. - for normalization, etc. 
    #     c0, c1, c2, c3, c4 = coefficientsList
    #     # maxHFTerm = 0 #c0 * self.fwd.filteredCalculateMaxHF(q_mesh_all, unconstrainedFaces)    #try not dividing by initial value
    #     # sumHFTerm = 0 #c1 * (self.fwd.calculateHFMeshSum(q_mesh_all) / numFaces) 
    #     # sumVertexDefects, maxVertexDefects, maxAngleBetweenNormals = calculateNormalsDiff(trimeshSolid)  
    #     imcTerm = c2 * self.calculateIntegralMeanCurvature(vertices, faces, face_adjacency, face_adjacency_edges)
    #     # imcTerm = c2 * (Solid.calculateSurfaceArea(vertices, faces) / initialParams[0])
    #     # normalsPenalty = c2 * sumVertexDefects
    #     vertexDefectsTerm = 0#c2 * calculateVertexDefects(vertices, faces, face_adjacency)
    #     maxNormalsTerm = 0#c3 * maxAngleBetweenNormals   
    #     #c4 was originally a thing but i've given up
    #     # return [vertexDefectsTerm + maxNormalsTerm, vertexDefectsTerm, maxNormalsTerm]
    #     #return [maxHFTerm + sumHFTerm + normalsPenalty + maxNormalsTerm + maxVertexDefectsTerm, normalsPenalty, maxNormalsTerm, maxVertexDefectsTerm]
    #     # return [normalsPenalty + maxNormalsTerm + maxVertexDefectsTerm, normalsPenalty, maxNormalsTerm]   
    #     return [imcTerm + maxNormalsTerm + vertexDefectsTerm, imcTerm, maxNormalsTerm]
    

    
