import numpy as np
import trimesh
from scipy.sparse import coo_matrix

class ObjectiveFunctionTools: 

    def __init__(self, chmod=0o774):
        #speed benchmarking tools
        self.testN = 0
        self.testTime = 0
        self.chmod = chmod
        return
    
    def setMeshAndGrids(self, trimeshSolid):
        self.all_faces = trimeshSolid.faces
        self.faces_sparse = trimeshSolid.faces_sparse
        self.face_adjacency = trimeshSolid.face_adjacency
        self.face_adjacency_edges = trimeshSolid.face_adjacency_edges
        vertices = trimeshSolid.vertices
        self.currentVerticesGrid = np.array([np.tile(vertices[np.newaxis, :], (len(vertices), 1, 1)) for _ in range(3)])
        self.newVerticesGrid = np.array([np.tile(vertices[np.newaxis, :], (len(vertices), 1, 1)) for _ in range(3)])
        self.initialIMC = self.calculateNonNormalizedIntegralMeanCurvature(vertices, self.all_faces, self.face_adjacency, self.face_adjacency_edges)
        self.initialVertexDefects = np.sum(np.abs(self.calculateVertexDefects(vertices, self.all_faces, self.face_adjacency)))
        # print(f"Shape of currentVerticesGrid: {self.currentVerticesGrid.shape}, shape of newVerticesGrid: {self.newVerticesGrid.shape}")
        # print(f"initial mesh triangles: {trimeshSolid.triangles}")
        # print(f"intial mesh triangle angle format: {trimeshSolid.triangles[100, 1]}")
        # input()
        # print(f"Vertex defects - trimesh: {trimeshSolid.vertex_defects}")
        print(f"Initial vertex defects sum: {self.initialVertexDefects}")
        return 
    
    def setParams(self, initialParams, coefficientsList):
        self.initialParams = initialParams
        self.coefficientsList = coefficientsList
        return 
    
    def setFacesToMove(self, facesToMove):
        self.facesToMove = facesToMove
        return 

    def setNewVerticesGrid(self, newVerticesGrid):
        self.newVerticesGrid = newVerticesGrid
        return 
    
    def calculateFaceNormals(self, vertices, faces):
        #get vectors of the two edges for each face
        vecs = vertices[faces[:, 1:]] - vertices[faces[:, 0, np.newaxis]]
        
        #calculate normal vectors using cross product
        faceNormals = np.cross(vecs[:, 0], vecs[:, 1])
        
        # Normalize each normal vector
        norms = np.linalg.norm(faceNormals, axis=1, keepdims=True)
        faceNormals /= norms
        
        return faceNormals
    
    ###OBJECTIVE FUNCTIONS AND CONSTRAINTS###

    def calculateNormalsDiff(self, trimeshSolid):
        # vertex_defects = trimesh.curvature.vertex_defects(trimeshSolid)
        vertex_defects = self.calculateVertexDefects(trimeshSolid.vertices, trimeshSolid.faces, trimeshSolid.face_adjacency)
        # vertex_defects = self.calculateVertexDefects(trimeshSolid.vertices, trimeshSolid.faces, trimeshSolid.face_adjacency)
        sumVertexDefects = np.sum(np.abs(vertex_defects))
        maxVertexDefect = np.max(np.abs(vertex_defects))  
        maxAngleBetweenNormals = 0 #TODO - will need to bring in these terms from impelmentation in RunModel
        return sumVertexDefects, maxVertexDefect, maxAngleBetweenNormals
    
    def calculateAngles(self, pairs):
        # pairs = np.asanyarray(pairs, dtype=np.float64)
        # do the dot product between vectors
        dots = np.sum(pairs[:, 0] * pairs[:, 1], axis=1)
        # clip for floating point error
        dots = np.clip(dots, -1.0, 1.0)
        # do cos and remove arbitrary sign
        angles = np.abs(np.arccos(dots))
        return angles

    def calculateDotProducts(self, vecA, vecB):
        #from trimesh - util.py
        vecA = np.asanyarray(vecA)
        # 3x faster than (a * b).sum(axis=1)
        # avoiding np.ones saves 5-10% sometimes
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
        # del vertices
        # del face_normals
        # del pairs 
        integralMeanCurvature = (angles * edges_length).sum() * 0.5
        return integralMeanCurvature
    
    def calculateIntegralMeanCurvature(self, vertices, faces, face_adjacency, face_adjacency_edges, initialIMC):
        """
        Compute the integral mean curvature of a mesh. Normalized: sum(curvatures / sum of curvatures)
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
        allValues = (angles * edges_length) * 0.5
        imc = np.sum(allValues)
        integralMeanCurvature = imc / initialIMC
        # del angles
        # del edges_length
        #integralMeanCurvature_all = (angles * edges_length).sum() * 0.5
        return integralMeanCurvature   
    
    def calculateIntegralMeanCurvatureParallel(self, vtx, dim):
        """
        Compute the integral mean curvature of a mesh. Normalized: sum(curvatures / sum of curvatures)
        """
        vertices = self.newVerticesGrid[dim][vtx]
        face_normals = self.calculateFaceNormals(vertices, self.all_faces)
        pairs = face_normals[self.face_adjacency]
        angles = self.calculateAngles(pairs)
        #maybe this could be added if we want to penalize concave angles more? but not sure if this could do more harm than good
        # angle_weights = np.where(angles > np.pi, 1.0, 0.5)
        # weighted_angles = angles * angle_weights
        edges_length = np.linalg.norm(np.subtract(
            *vertices[self.face_adjacency_edges.T]), axis=1)
        integralMeanCurvature = np.sum((angles * edges_length) * 0.5) / self.initialIMC
        # edges_length = np.linalg.norm(np.subtract(
        #     *vertices[self.face_adjacency_edges.T]), axis=1)
        # integralMeanCurvature = np.sum((weighted_angles * edges_length) * 0.5) / self.initialIMC
        return integralMeanCurvature   
    

    def calcUnitVectors(self, vectors, threshold=None):
        #from trimesh - utils - unitize
        vectors = np.asanyarray(vectors)
        # allow user to set zero threshold
        if threshold is None:
            threshold = np.finfo(np.float64).resolution * 100
        if len(vectors.shape) == 2:
            # for (m, d) arrays take the per-row unit vector
            # using sqrt and avoiding exponents is slightly faster
            # also dot with ones is faser than .sum(axis=1)
            norm = np.sqrt(np.dot(vectors * vectors,
                                [1.0] * vectors.shape[1]))
            # non-zero norms
            valid = norm > threshold
            # in-place reciprocal of nonzero norms
            norm[valid] **= -1
            # multiply by reciprocal of norm
            unitVectors = vectors * norm.reshape((-1, 1))

        elif len(vectors.shape) == 1:
            # treat 1D arrays as a single vector
            norm = np.sqrt(np.dot(vectors, vectors))
            valid = norm > threshold
            if valid:
                unitVectors = vectors / norm
            else:
                unitVectors = vectors.copy()
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

        anglesMatrixSparse = coo_matrix((
                angles.flatten(),
                (self.faces_sparse.row, self.faces_sparse.col)),
                self.faces_sparse.shape)

        # angle_sum = np.array(angles.sum(axis=1)).flatten()
        angle_sum = np.array(anglesMatrixSparse.sum(axis=1)).flatten()
        # print(f"defects: {(2*np.pi) - angle_sum}")
        # input()

        return (2*np.pi) - angle_sum
    
    def calculateVertexDefectsParallel(self, vtx, dim):
        #calculate mesh face angles - so angles of triangles
        vertices = self.newVerticesGrid[dim][vtx]
        triangles = np.asanyarray(vertices.view(np.ndarray)[self.all_faces], dtype=np.float64)
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

        anglesMatrixSparse = coo_matrix((
                angles.flatten(),
                (self.faces_sparse.row, self.faces_sparse.col)),
                self.faces_sparse.shape)

        # angle_sum = np.array(angles.sum(axis=1)).flatten()
        angle_sum = np.array(anglesMatrixSparse.sum(axis=1)).flatten()
        # print(f"defects: {(2*np.pi) - angle_sum}")
        return (2*np.pi) - angle_sum
    

# angle_sum = np.array(mesh.face_angles_sparse.sum(axis=1)).flatten()
# defect = (2 * np.pi) - angle_sum
# return defect

#mesh.face_angles_sparse: flattened then sparsed mesh.face_angles
#mesh.face_angles: triangles.angles(self.triangles)
#triangles.angles: degenerate angles will be 0
# # don't copy triangles
# triangles = np.asanyarray(triangles, dtype=np.float64)

# # get a unit vector for each edge of the triangle
# u = unitize(triangles[:, 1] - triangles[:, 0])
# v = unitize(triangles[:, 2] - triangles[:, 0])
# w = unitize(triangles[:, 2] - triangles[:, 1])

# # run the cosine and per-row dot product
# result = np.zeros((len(triangles), 3), dtype=np.float64)
# # clip to make sure we don't float error past 1.0
# result[:, 0] = np.arccos(np.clip(diagonal_dot(u, v), -1, 1))
# result[:, 1] = np.arccos(np.clip(diagonal_dot(-u, w), -1, 1))
# # the third angle is just the remaining
# result[:, 2] = np.pi - result[:, 0] - result[:, 1]

# # a triangle with any zero angles is degenerate
# # so set all of the angles to zero in that case
# result[(result < tol.merge).any(axis=1), :] = 0.0

#unitize: 
#tol.merge: 
#diagonal_dot: 



    # def vtxFacesObjectiveFunctionCalc(self, verticesList):
    #     c0, c1, c2, c3, c4 = self.coefficientsList
    #     imcTerm = c2 * self.calculateIntegralMeanCurvature(verticesList, self.all_faces, self.face_adjacency, self.face_adjacency_edges, self.initialIMC)
    #     vertexDefectsTerm = 0#c2 * calculateVertexDefects(vertices, faces, face_adjacency)
    #     maxNormalsTerm = 0#c3 * maxAngleBetweenNormals   
    #     return imcTerm + maxNormalsTerm + vertexDefectsTerm #[imcTerm + maxNormalsTerm + vertexDefectsTerm, imcTerm, maxNormalsTerm]

    def vtxFacesObjectiveFunctionCalc(self, verticesList):
        c0, c1, c2, c3, c4 = self.coefficientsList
        imcTerm = c2 * self.calculateIntegralMeanCurvature(verticesList, self.all_faces, self.face_adjacency, self.face_adjacency_edges, self.initialIMC)
        vertexDefectsTerm = c3 * (np.sum(np.abs(self.calculateVertexDefects(verticesList, self.all_faces, self.face_adjacency))) / self.initialVertexDefects)
        maxNormalsTerm = 0#c3 * maxAngleBetweenNormals   
        return imcTerm + maxNormalsTerm + vertexDefectsTerm #[imcTerm + maxNormalsTerm + vertexDefectsTerm, imcTerm, maxNormalsTerm]

    # def objectiveFunction(self, vtx, dim): 
    #     c0, c1, c2, c3, c4 = self.coefficientsList
    #     imcTerm = c2 * self.calculateIntegralMeanCurvatureParallel(vtx, dim)
    #     vertexDefectsTerm = 0#c2 * calculateVertexDefects(vertices, faces, face_adjacency)
    #     maxNormalsTerm = 0#c3 * maxAngleBetweenNormals   
    #     return imcTerm + maxNormalsTerm + vertexDefectsTerm #[imcTerm + maxNormalsTerm + vertexDefectsTerm, imcTerm, maxNormalsTerm]

    def objectiveFunction(self, vtx, dim): 
        c0, c1, c2, c3, c4 = self.coefficientsList
        imcTerm = c2 * self.calculateIntegralMeanCurvatureParallel(vtx, dim)
        vertexDefectsTerm = c3 * (np.sum(np.abs(self.calculateVertexDefectsParallel(vtx, dim))) / self.initialVertexDefects)
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
    

    
