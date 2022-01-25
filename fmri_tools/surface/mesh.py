# -*- coding: utf-8 -*-

import functools
import numpy as np
from scipy.sparse import csr_matrix, triu, dia_matrix

__all__ = ['Mesh']


class Mesh:
    """Implementation of some useful functions to work with a triangle surface
    mesh.

    Parameters
    ----------
    verts : np.ndarray, shape=(N,3)
        Vertex coordinates.
    faces : np.ndarray, shape=(M,3)
        Vertex indices of each triangle.

    Raises
    ------
    ValueError :
        If `verts` has a wrong shape or `faces` does not match the vertex array.

    """

    def __init__(self, verts, faces):
        self.verts = verts
        self.faces = faces

    @property
    @functools.lru_cache
    def tris(self):
        """Array of coordinates in each face.

        Returns
        -------
        np.ndarray, shape=(N,3,3)
            Three coordinates for each face in the mesh which consists of three
            points.

        """

        return self.verts[self.faces]

    @property
    @functools.lru_cache
    def adjm(self):
        """Compute a sparse adjacency matrix. The matrix has the size
        (nvertex, nvertex). Each matrix entry with value 1 stands for an edge of
        the surface mesh.

        Returns
        -------
        scipy.sparse.csr.csr_matrix
            Sparse adjacency matrix.

        """

        # number of vertices
        nverts = len(self.verts)

        # initialise
        row = []
        col = []

        # get rows and columns of edges
        row.extend(list(self.faces[:, 0]))
        col.extend(list(self.faces[:, 1]))

        row.extend(list(self.faces[:, 1]))
        col.extend(list(self.faces[:, 2]))

        row.extend(list(self.faces[:, 2]))
        col.extend(list(self.faces[:, 0]))

        # make sure that all edges are symmetric
        row.extend(list(self.faces[:, 1]))
        col.extend(list(self.faces[:, 0]))

        row.extend(list(self.faces[:, 2]))
        col.extend(list(self.faces[:, 1]))

        row.extend(list(self.faces[:, 0]))
        col.extend(list(self.faces[:, 2]))

        # adjacency entries get value 1
        data = np.ones(len(row), dtype=np.int8)

        return csr_matrix((data, (row, col)), shape=(nverts, nverts))

    @property
    @functools.lru_cache
    def edges(self):
        """Edge array.

        From the adjacency matrix, all edges in the surface mesh are extracted.
        Only the upper triangular portion of the matrix is used to avoid double
        countings of edges.

        Returns
        -------
        edge_coords : np.ndarray, shape=(N,2)
            Vertex indices of each edge.

        """

        adjm_upper = triu(self.adjm, 1)
        edge_coords = np.zeros((len(adjm_upper.row), 2), dtype=np.int64)
        edge_coords[:, 0] = adjm_upper.row
        edge_coords[:, 1] = adjm_upper.col

        return edge_coords

    @property
    @functools.lru_cache
    def vfm(self):
        """Compute a sparse matrix of vertex-face associations. The matrix has
        the size (nvertex, nfaces). For each vertex index, all associated faces
        are listed.

        Returns
        -------
        scipy.sparse.csr.csr_matrix
            Sparse adjacency matrix.

        """

        # number of vertices and faces
        nverts = len(self.verts)
        nfaces = len(self.faces)

        row = np.hstack(self.faces.T)
        col = np.tile(range(nfaces), (1, 3)).squeeze()

        # vertex-face associations get value 1
        data = np.ones(len(row), dtype=np.int8)

        return csr_matrix((data, (row, col)), shape=(nverts, nfaces))

    @property
    @functools.lru_cache
    def face_normals(self):
        """Face-wise surfaces normals.

        Returns
        -------
        n : np.ndarray, shape=(M,3)
            Array of face-wise normal vectors.

        """

        # calculate the normal for all triangles by taking the cross product of
        # the vectors v1-v0 and v2-v0 in each triangle and normalize
        n = np.cross(self.tris[::, 1] - self.tris[::, 0],
                     self.tris[::, 2] - self.tris[::, 0])
        n = self._normalize_array(n)

        return n

    @property
    @functools.lru_cache
    def vertex_normals(self):
        """Vertex-wise surfaces normals. The code is taken from [1]_ and adapted
        to my own purposes.

        Returns
        -------
        norm : np.ndarray, shape=(N,)
            Array of vertex-wise normal vectors.

        References
        -------
        .. [1] https://sites.google.com/site/dlampetest/python/calculating-normals-of-a-triangle-mesh-using-numpy

        """

        # face normals
        n = self.face_normals

        # calculate vertex-wise normals from face normals and normalize
        n = self._f2v(n)
        n = self._normalize_array(n)

        return n

    @property
    @functools.lru_cache
    def face_areas(self):
        """Triangle areas.

        Returns
        -------
        n : np.ndarray, shape=(M,)
            Array of face areas.

        """

        # calculate the normal for all triangles by taking the cross product of
        # the vectors v1-v0 and v2-v0 in each triangle and get face area from
        # length
        n = np.cross(self.tris[::, 1] - self.tris[::, 0],
                     self.tris[::, 2] - self.tris[::, 0])
        n = np.sqrt((n ** 2).sum(-1)) / 2

        return n

    @property
    @functools.lru_cache
    def cotangent(self):
        """Cotangent of angle opposite each vertex in each face.

        Returns
        -------
        cots : np.ndarray, shape=(3,M)
            Cotangent for each vertex in each face.

        """

        tris10 = self.tris[:, 1] - self.tris[:, 0]
        tris20 = self.tris[:, 2] - self.tris[:, 0]
        cots0 = (tris10 * tris20).sum(1) / \
            np.sqrt((np.cross(tris10, tris20) ** 2).sum(1))

        tris21 = self.tris[:, 2] - self.tris[:, 1]
        tris01 = self.tris[:, 0] - self.tris[:, 1]
        cots1 = (tris21 * tris01).sum(1) / \
            np.sqrt((np.cross(tris21, tris01) ** 2).sum(1))

        tris02 = self.tris[:, 0] - self.tris[:, 2]
        tris12 = self.tris[:, 1] - self.tris[:, 2]
        cots2 = (tris02 * tris12).sum(1) / \
            np.sqrt((np.cross(tris02, tris12) ** 2).sum(1))

        # remove invalid values
        cots = np.vstack([cots0, cots1, cots2])
        cots[np.isinf(cots)] = 0
        cots[np.isnan(cots)] = 0

        return cots

    @property
    @functools.lru_cache
    def laplace_beltrami(self):
        """Laplace-Beltrami operator. The operator is a sparse adjacency matrix
        with edge weights determined by the cotangents of the angles opposite
        each edge. See [1]_ and [2]_ for further details.

        Returns
        -------
        scipy.sparse.csr.csr_matrix
            Laplace-Beltrami operator.

        References
        ----------
        .. [1] Chen, Yi, et al. Scale-specific analysis of fMRI data on the
        irregular cortical surface. Neuroimage 181, 370--381 (2018).
        .. [2] https://gallantlab.github.io/generated/cortex.polyutils.Surface.html

        """

        nverts = len(self.verts)
        cots0, cots1, cots2 = self.cotangent

        # D: Diagonal elements of the lumped mass matrix. For each vertex, the
        # sum over areas of all assocated face areas is computed.
        D_diag = self.vfm.dot(self.face_areas) / 3.0
        D_inv = dia_matrix((1.0 / D_diag, [0]), (nverts, nverts))

        # W is weighted adjacency matrix
        W0 = csr_matrix(
            (cots0, (self.faces[:, 1], self.faces[:, 2])), (nverts, nverts))
        W1 = csr_matrix(
            (cots1, (self.faces[:, 2], self.faces[:, 0])), (nverts, nverts))
        W2 = csr_matrix(
            (cots2, (self.faces[:, 0], self.faces[:, 1])), (nverts, nverts))
        W = (W0 + W0.T + W1 + W1.T + W2 + W2.T) / 2.0

        # V is a diagonal matrix that normalizes the adjacencies
        V = dia_matrix((np.array(W.sum(0)).ravel(), [0]), (nverts, nverts))

        return D_inv * (V - W)

    @property
    @functools.lru_cache
    def avg_edge_length(self):
        """Average of all edges in the surface.
        
        Returns
        -------
        float
            Average edge length.

        """
        
        edge_length = np.sqrt(((self.verts[self.edges[:, 0]] -
                               self.verts[self.edges[:, 1]])**2).sum(axis=1))

        return edge_length.mean()

    @property
    @functools.lru_cache
    def n_neighbors(self):
        """Number of neighbors for each vertex.

        The returned array is divided by max() to account for the fact that the
        entries in the sparse matrix are not necessarily ones.

        Returns
        -------
        np.ndarray, shape=(N,)
            Array of vertex neighbors.

        """

        return self.adjm.sum(axis=0) / self.adjm.max()

    def neighborhood(self, ind):
        """Compute 1-ring neighborhood for one vertex.

        Parameters
        ----------
        ind : int
            Vertex index.

        Returns
        -------
        np.ndarray, shape=(N,)
            Array of neighborhood indices.

        """

        return self.adjm[ind, :].indices

    def _f2v(self, nf_arr):
        """Transform face- to vertex-wise expression.

        Parameters
        ----------
        nf_arr : np.ndarray, shape=(M,3)
            Face-wise array.

        Returns
        -------
        np.ndarray, shape=(N,3)
            Vertex-wise array.

        """

        return self.vfm.dot(nf_arr)

    @staticmethod
    def _normalize_array(arr):
        """Normalize a numpy array of shape=(n,3) along axis=1.

        Parameters
        ----------
        arr : np.ndarray, shape=(N,3)
            Data array

        Returns
        -------
        res : np.ndarray, shape=(N,3)
            Normalized data array.

        """

        lens = np.sqrt(arr[:, 0] ** 2 + arr[:, 1] ** 2 + arr[:, 2] ** 2)
        lens[lens == 0] = np.nan
        res = np.zeros_like(arr)
        res[:, 0] = arr[:, 0] / lens
        res[:, 1] = arr[:, 1] / lens
        res[:, 2] = arr[:, 2] / lens
        res_sum = np.sum(res, axis=1)
        res[~np.isfinite(res_sum), :] = 0

        return res

    @property
    def verts(self):
        return self._verts

    @verts.setter
    def verts(self, v):
        v = np.asarray(v)
        if v.ndim != 2 or np.shape(v)[1] != 3:
            raise ValueError("Vertices have wrong shape!")

        self._verts = v

    @property
    def faces(self):
        return self._faces

    @faces.setter
    def faces(self, f):
        f = np.asarray(f)
        if f.ndim != 2 or np.shape(f)[1] != 3:
            raise ValueError("Vertices have wrong shape!")

        if np.max(f) != len(self.verts) - 1:
            raise ValueError("Faces do not match vertex array!")

        self._faces = f
