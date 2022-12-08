# -*- coding: utf-8 -*-

import functools
import numpy as np
import nibabel as nb
from nibabel.freesurfer.io import read_geometry
from scipy.sparse import csr_matrix, triu, dia_matrix

# local inputs
from ..io.affine import read_vox2ras_tkr
from ..utils.interpolation import linear_interpolation3d
from ..utils.apply_affine_chunked import apply_affine_chunked

__all__ = ["Mesh"]


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
    @functools.lru_cache()
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
    @functools.lru_cache()
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
    @functools.lru_cache()
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
    @functools.lru_cache()
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
    @functools.lru_cache()
    def face_normals(self):
        """Face-wise surfaces normals.

        Returns
        -------
        n : np.ndarray, shape=(M,3)
            Array of face-wise normal vectors.

        """

        # calculate the normal for all triangles by taking the cross product of
        # the vectors v1-v0 and v2-v0 in each triangle and normalize
        n = np.cross(
            self.tris[::, 1] - self.tris[::, 0], self.tris[::, 2] - self.tris[::, 0]
        )
        n = self._normalize_array(n)

        return n

    @property
    @functools.lru_cache()
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
    @functools.lru_cache()
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
        n = np.cross(
            self.tris[::, 1] - self.tris[::, 0], self.tris[::, 2] - self.tris[::, 0]
        )
        n = np.sqrt((n**2).sum(-1)) / 2

        return n

    @property
    @functools.lru_cache()
    def cotangent(self):
        """Cotangent of angle opposite each vertex in each face.

        Returns
        -------
        cots : np.ndarray, shape=(3,M)
            Cotangent for each vertex in each face.

        """

        tris10 = self.tris[:, 1] - self.tris[:, 0]
        tris20 = self.tris[:, 2] - self.tris[:, 0]
        cots0 = (tris10 * tris20).sum(1) / np.sqrt(
            (np.cross(tris10, tris20) ** 2).sum(1)
        )

        tris21 = self.tris[:, 2] - self.tris[:, 1]
        tris01 = self.tris[:, 0] - self.tris[:, 1]
        cots1 = (tris21 * tris01).sum(1) / np.sqrt(
            (np.cross(tris21, tris01) ** 2).sum(1)
        )

        tris02 = self.tris[:, 0] - self.tris[:, 2]
        tris12 = self.tris[:, 1] - self.tris[:, 2]
        cots2 = (tris02 * tris12).sum(1) / np.sqrt(
            (np.cross(tris02, tris12) ** 2).sum(1)
        )

        # remove invalid values
        cots = np.vstack([cots0, cots1, cots2])
        cots[np.isinf(cots)] = 0
        cots[np.isnan(cots)] = 0

        return cots

    @property
    @functools.lru_cache()
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
        .. [2] https://gallantlab.github.io/pycortex/generated/cortex.polyutils.Surface.html

        """

        nverts = len(self.verts)
        cots0, cots1, cots2 = self.cotangent

        # D: Diagonal elements of the lumped mass matrix. For each vertex, the
        # sum over areas of all assocated face areas is computed.
        D_diag = self.vfm.dot(self.face_areas) / 3.0
        D_inv = dia_matrix((1.0 / D_diag, [0]), (nverts, nverts))

        # W is weighted adjacency matrix
        W0 = csr_matrix((cots0, (self.faces[:, 1], self.faces[:, 2])), (nverts, nverts))
        W1 = csr_matrix((cots1, (self.faces[:, 2], self.faces[:, 0])), (nverts, nverts))
        W2 = csr_matrix((cots2, (self.faces[:, 0], self.faces[:, 1])), (nverts, nverts))
        W = (W0 + W0.T + W1 + W1.T + W2 + W2.T) / 2.0

        # V is a diagonal matrix that normalizes the adjacencies
        V = dia_matrix((np.array(W.sum(0)).ravel(), [0]), (nverts, nverts))

        return D_inv * (V - W)

    @property
    @functools.lru_cache()
    def boundary_vertices(self):
        """Determination of boundary vertices in surface mesh. The implementation
        follows [1]_. The algorithm uses the property that every edge appears in either
        one (border edge) or two (non-border edge) faces.

        Returns
        -------
        ndarray, shape=(N,)
            Array containing vertex indices of border vertices.

        References
        ----------
        .. [1] https://gallantlab.github.io/pycortex/generated/cortex.polyutils.Surface.html

        """

        first = np.hstack(
            [
                self.faces[:, 0],
                self.faces[:, 1],
                self.faces[:, 2],
            ]
        )
        second = np.hstack(
            [
                self.faces[:, 1],
                self.faces[:, 2],
                self.faces[:, 0],
            ]
        )

        polygon_edges = np.vstack([first, second])
        polygon_edges = np.vstack(
            [polygon_edges.min(axis=0), polygon_edges.max(axis=0)]
        )

        sort_order = np.lexsort(polygon_edges)
        sorted_edges = polygon_edges[:, sort_order]

        duplicate_mask = (sorted_edges[:, :-1] == sorted_edges[:, 1:]).sum(axis=0) == 2

        nonduplicate_mask = np.ones(sorted_edges.shape[1], dtype=bool)
        nonduplicate_mask[:-1][duplicate_mask] = False
        nonduplicate_mask[1:][duplicate_mask] = False

        border_label = np.hstack(
            [
                sorted_edges[:, nonduplicate_mask][0, :],
                sorted_edges[:, nonduplicate_mask][1, :],
            ]
        )
        border_label = np.unique(border_label)

        return border_label

    @property
    @functools.lru_cache()
    def avg_edge_length(self):
        """Average of all edges in the surface.

        Returns
        -------
        float
            Average edge length.

        """

        edge_length = np.sqrt(
            ((self.verts[self.edges[:, 0]] - self.verts[self.edges[:, 1]]) ** 2).sum(
                axis=1
            )
        )

        return edge_length.mean()

    @property
    @functools.lru_cache()
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

    def transform_coords(self, file_cmap, file_target):
        """Transform vertex coordinates to a target space using a coordinate mapping
        file. Vertices outside the coordinate mapping are set to nan.

        Parameters
        ----------
        file_cmap : str
            File name of coordinate mapping.
        file_target : str
            File name of volume in target space.

        Returns
        -------
        np.ndarray, shape=(N,3)
            Vertex-wise array.

        """

        nx, ny, nz = nb.load(file_cmap).header["dim"][1:4]
        _, ras2vox = read_vox2ras_tkr(file_cmap)
        vox2ras, _ = read_vox2ras_tkr(file_target)

        # transform to voxel space
        vtx_vox = apply_affine_chunked(ras2vox, self.verts)

        # mask voxels which are outside the coordinate map
        mask = np.ones(len(self.verts), dtype=bool)
        mask[vtx_vox[:, 0] < 0] = 0
        mask[vtx_vox[:, 1] < 0] = 0
        mask[vtx_vox[:, 2] < 0] = 0
        mask[vtx_vox[:, 0] > nx - 1] = 0
        mask[vtx_vox[:, 1] > ny - 1] = 0
        mask[vtx_vox[:, 2] > nz - 1] = 0
        mask = mask[mask == 1]
        vtx_vox = vtx_vox[mask, :]

        # apply transformation
        arr_cmap = nb.load(file_cmap).get_fdata()
        x = linear_interpolation3d(vtx_vox[:, 0], vtx_vox[:, 1], vtx_vox[:, 2],
                                   arr_cmap[:, :, :, 0])
        y = linear_interpolation3d(vtx_vox[:, 0], vtx_vox[:, 1], vtx_vox[:, 2],
                                   arr_cmap[:, :, :, 1])
        z = linear_interpolation3d(vtx_vox[:, 0], vtx_vox[:, 1], vtx_vox[:, 2],
                                   arr_cmap[:, :, :, 2])

        # update vertex array
        vtx_res = apply_affine_chunked(vox2ras, np.array([x, y, z]).T)
        self.verts = np.empty_like(self.verts)
        self.verts[:] = np.nan
        self.verts[mask, :] = vtx_res

        return self.verts

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

    @classmethod
    def from_file(cls, file_surf):
        """Initialize class object from file.

        Parameters
        ----------
        file_surf : str
            File name of freesurfer geometry.

        Returns
        -------
        Constructor from file name.

        """

        vtx, fac = read_geometry(file_surf)

        return cls(vtx, fac)

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
