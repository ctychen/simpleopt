import numpy as np
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
        q_mesh_initial = self.calculateAllHF_AllVals_GivenVertices(self.q_dir, self.q_mag, vertices)
        self.maxHF_initial = self.filteredCalculateMaxHF(q_mesh_initial, unconstrainedFaces=[]) #self.fwd.filteredCalculateMaxHF(q_mesh_initial, unconstrainedFaces = [])
        self.sumHF_initial = self.calculateHFMeshSum(q_mesh_initial) 
        self.initialHFSum = self.calculateHFMeshSum(q_mesh_initial)
        self.initialVertexDefects = np.sum(np.abs(self.calculateVertexDefects(vertices, self.all_faces, self.face_adjacency)))
        print(f"Initial max HF: {self.maxHF_initial}")
        print(f"Initial sum HF: {self.sumHF_initial}")
        print(f"Initial vertex defects sum: {self.initialVertexDefects}")
        print(f"Initial integral mean curvature: {self.initialIMC}")
        return 
    
    def setForwardModel(self, fwdModel):
        self.fwdModel = fwdModel
        return
    
    def setParams(self, initialParams, coefficientsList):
        self.initialParams = initialParams
        self.coefficientsList = coefficientsList
        return 
    
    def setHeatFluxParams(self, q_dir, q_mag):
        self.q_dir = q_dir
        self.q_mag = q_mag
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

    def calculateNormalsDiff(self, trimeshSolid):
        vertex_defects = self.calculateVertexDefects(trimeshSolid.vertices, trimeshSolid.faces, trimeshSolid.face_adjacency)
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
        return integralMeanCurvature   
    

    def calcUnitVectors(self, vectors, threshold=None):
        #reference from trimesh - utils - unitize
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

    def calculateVertexDefects(self, vertices, faces):
        #reference from Trimesh - curvature - vertex defects
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

        angle_sum = np.array(anglesMatrixSparse.sum(axis=1)).flatten()

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

        angle_sum = np.array(anglesMatrixSparse.sum(axis=1)).flatten()
        return (2*np.pi) - angle_sum
    

    def calculateAllHF_AllVals(self, q_dir, q_mag, vtx, dim):
        """
        Calculate HF on every mesh element, with no exclusions, and returning them all in a list
        """
        vertices = self.newVerticesGrid[dim][vtx]
        normals = self.calculateFaceNormals(vertices, self.all_faces)
        q_mesh_all = -1 * (np.dot(normals, q_dir)) * q_mag
        return q_mesh_all    
    
    def calculateAllHF_AllVals_GivenVertices(self, q_dir, q_mag, vertices):
        """
        Calculate HF on every mesh element, with no exclusions, and returning them all in a list
        """
        normals = self.calculateFaceNormals(vertices, self.all_faces)
        q_mesh_all = -1 * (np.dot(normals, q_dir)) * q_mag
        return q_mesh_all    
    
    def calculateAllHF(self, hfMode, q_dir, q_mag, vtx, dim):
        """
        calculate HF on every mesh element, but this should use nonuniform HF on mesh element centers
        this could also use uniform HF; just need meshHFProfile to be a list of the same HF value for each element
        """ 

        q_mesh_all = []

        if hfMode == 'uniform':
            q_mesh_all = self.calculateAllHF_AllVals(q_dir, q_mag, vtx, dim)

        elif hfMode == "uniform_multiple":
            #this works for 2 HF magnitudes and directions but not more - would have to generalize?
            vertices = self.newVerticesGrid[dim][vtx]
            normals = self.calculateFaceNormals(vertices, self.all_faces)
            dot_product_q0 = np.dot(normals, q_dir[0])
            mask_q0 = dot_product_q0 > 0
            q0_mesh = -1 * np.where(mask_q0, 0, dot_product_q0) * q_mag[0]

            dot_product_q1 = np.dot(normals, q_dir[1])
            mask_q1 = dot_product_q1 > 0
            q1_mesh = -1 * np.where(mask_q1, 0, dot_product_q1) * q_mag[1]

            q_mesh_all = np.add(q0_mesh, q1_mesh)

        return q_mesh_all
    
    def calculateHFMeshSum(self, q_mesh_all):
        """
        Calculate sum of heat flux from all mesh elements
        """
        q_mesh_vals = np.where(q_mesh_all > 0, q_mesh_all, 0)
        q_mesh_sum = np.sum(np.abs(q_mesh_vals))
        return q_mesh_sum
    
    def calculateMaxHF(self, q_mesh_all):
        return np.max(q_mesh_all)

    def filteredCalculateMaxHF(self, q_mesh_all, unconstrainedFaces = []):
        """
        filtering out unconstrained faces from max HF calculation - originally was added bc of HF with incident angle
        if no unconstrained faces defined, then we just return max HF same as before
        """

        # if unconstrainedFaces: #if unconstrainedFaces is not empty
        #     unconstrainedFaces = list(unconstrainedFaces) 
        #     mask = np.ones(q_mesh_all.shape, dtype=bool)
        #     mask[unconstrainedFaces] = False
        #     q_mesh_all[mask] = 0

        return np.max(q_mesh_all)

    def calculateHFDistribution(self, q_dir, q_mag, vtx, dim):
        """
        Calculate mean, variance, standard deviation on all HF over mesh elements
        """
        vertices = self.newVerticesGrid[dim][vtx]
        normals = self.calculateFaceNormals(vertices, self.all_faces)
        q_mesh_all = -1 * (np.dot(normals, q_dir)) * q_mag

        hfMean = np.mean(q_mesh_all)
        hfVariance = np.var(q_mesh_all)
        hfStd = np.std(q_mesh_all)

        return hfMean, hfVariance, hfStd, q_mesh_all

    def vtxFacesObjectiveFunctionCalc(self, verticesList):
        c0, c1, c2, c3, c4 = self.coefficientsList
        #self.hfMode, self.q_dir, self.q_mag, self.facesToMove, 0), self.facesToMove
        q_mesh_all = self.calculateAllHF_AllVals_GivenVertices(self.fwdModel.q_dir, self.fwdModel.q_mag, verticesList)
        maxHeatFluxTerm = c0 * self.filteredCalculateMaxHF(q_mesh_all)
        sumHeatFluxTerm = c1 * (self.calculateHFMeshSum(q_mesh_all) / self.initialHFSum)
        imcTerm = c2 * self.calculateIntegralMeanCurvature(verticesList, self.all_faces, self.face_adjacency, self.face_adjacency_edges, self.initialIMC)
        vertexDefectsTerm = c3 * np.sum(np.abs(self.calculateVertexDefects(verticesList, self.all_faces)))
        return maxHeatFluxTerm + sumHeatFluxTerm + imcTerm + vertexDefectsTerm 

    def objectiveFunction(self, vtx, dim): 
        c0, c1, c2, c3, c4 = self.coefficientsList
        q_mesh_all = self.calculateAllHF(self.fwdModel.hfMode, self.fwdModel.q_dir, self.fwdModel.q_mag, self.facesToMove, vtx, dim)
        maxHeatFluxTerm = c0 * self.filteredCalculateMaxHF(q_mesh_all)
        sumHeatFluxTerm = c1 * (self.calculateHFMeshSum(q_mesh_all) / self.initialHFSum)
        imcTerm = c2 * self.calculateIntegralMeanCurvatureParallel(vtx, dim)
        vertexDefectsTerm = c3 * np.sum(np.abs(self.calculateVertexDefectsParallel(vtx, dim)))
        return maxHeatFluxTerm + sumHeatFluxTerm + imcTerm + vertexDefectsTerm
    
    
