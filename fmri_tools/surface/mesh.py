# -*- coding: utf-8 -*-
"""Mesh utilities."""

import functools
import itertools
import os

import nibabel as nb
import numpy as np
from nibabel.freesurfer.io import read_geometry, read_label, write_morph_data
from numpy.linalg import norm
from scipy.sparse import csr_matrix, dia_matrix, triu

from ..io.affine import read_vox2ras_tkr
from ..io.filename import get_filename
from ..io.surf import read_mgh, write_mgh
from ..registration.transform import apply_affine_chunked
from ..utils.interpolation import linear_interpolation3d

__all__ = [
    "Mesh",
    "calculate_area",
    "gradient",
    "b0_orientation",
    "clip_surface",
    "clip_mgh",
]


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
        for i, n in enumerate((nx, ny, nz)):
            mask[vtx_vox[:, i] < 0] = 0
            mask[vtx_vox[:, i] > n - 1] = 0
        vtx_vox = vtx_vox[mask == 1, :]

        # apply transformation
        arr_cmap = nb.load(file_cmap).get_fdata()
        x = linear_interpolation3d(
            vtx_vox[:, 0], vtx_vox[:, 1], vtx_vox[:, 2], arr_cmap[:, :, :, 0]
        )
        y = linear_interpolation3d(
            vtx_vox[:, 0], vtx_vox[:, 1], vtx_vox[:, 2], arr_cmap[:, :, :, 1]
        )
        z = linear_interpolation3d(
            vtx_vox[:, 0], vtx_vox[:, 1], vtx_vox[:, 2], arr_cmap[:, :, :, 2]
        )

        # update vertex array
        vtx_res = apply_affine_chunked(vox2ras, np.array([x, y, z]).T)
        self.verts = np.empty_like(self.verts)
        self.verts[:] = np.nan
        self.verts[mask, :] = vtx_res

        return self.verts

    def remove_vertices(self, ind_keep, create_ind=False):
        """Remove vertices and update the corresponding face array.

        Parameters
        ----------
        ind_keep : list
            Index list of vertices to keep.
        create_ind : bool
            Return cleaned index list.

        Returns
        -------
        np.ndarray, shape=(N,3)
            Array of remaining vertices.
        np.ndarray, shape=(M,3)
            Array of updated faces.
        list
            Returned only if `create_ind` is True. Updated index list of vertices to
            keep after vertex cleaning.

        """
        # get indices which will be removed
        ind_tmp = np.arange(len(self.verts))
        ind_remove = list(set(ind_tmp) - set(ind_keep))
        ind_remove = sorted(ind_remove, reverse=True)

        # get new vertices
        vtx = self.verts[ind_keep, :]

        # get new faces
        fac_keep = np.zeros(len(self.faces))
        fac_keep += np.in1d(self.faces[:, 0], ind_keep)
        fac_keep += np.in1d(self.faces[:, 1], ind_keep)
        fac_keep += np.in1d(self.faces[:, 2], ind_keep)
        fac = self.faces[fac_keep == 3, :]

        # reindex faces
        loop_status = 0
        loop_length = len(ind_remove)
        for i in range(loop_length):
            # print status
            counter = np.floor(i / loop_length * 100)
            if counter != loop_status:
                print("sort faces: " + str(counter) + " %")
                loop_status = counter

            tmp = fac[fac >= ind_remove[i]] - 1
            fac[fac >= ind_remove[i]] = tmp

        # get indices which will be cleaned
        ind_vtx = np.arange(len(vtx))
        ind_fac = list(itertools.chain(*fac))
        ind_fac = list(set(ind_fac))
        ind_remove = list(set(ind_vtx) - set(ind_fac))
        ind_remove = sorted(ind_remove, reverse=True)

        # remove singularities (vertices without faces)
        loop_status = 0
        loop_length = len(ind_remove)
        for i in range(loop_length):
            # print status
            counter = np.floor(i / loop_length * 100)
            if counter != loop_status:
                print("clean faces: " + str(counter) + " %")
                loop_status = counter

            # remove vertex and index
            vtx = np.delete(vtx, ind_remove[i], 0)
            ind_keep = np.delete(ind_keep, ind_remove[i], 0)

            # sort faces
            tmp = fac[fac >= ind_remove[i]] - 1
            fac[fac >= ind_remove[i]] = tmp

        self.verts = np.array(vtx)
        self.faces = np.array(fac)
        if not create_ind:
            return self.verts, self.faces
        return self.verts, self.faces, ind_keep

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


def calculate_area(filename_surf, filename_area=""):
    """The function calculates vertex-wise surface area. The code is taken from the
    octave code surf2area.m from Anderson Winkler found in his github repository [1].
    Consider a triangular face ABC with corner points
    a = [x_A, y_A, z_A]'
    b = [x_B, y_B, z_B]'
    c = [x_C, y_C, z_C]'
    The area for this triangle is given by the normed cross product A = |u x v|/2 with
    u = a - c and v = b - c. This is a face-wise surface area representation. To convert
    this to a vertex-wise representation, we assign each vertex one third of the sum of
    the areas of all faces that meet at that vertex, cf. [2].

    Parameters
    ----------
    filename_surf : str
        Input file geometry on which surface area is calculated.
    filename_area : str, optional
        File name of the surface area file. The default is "".

    Returns
    -------
    dpv : ndarray
        Vertex-wise surface area.

    References
    -------
    .. [1] https://github.com/andersonwinkler/areal
    .. [2] Winkler, A, et al. Measuring and comparing brain cortical surface
    area and other areal quantities, Neuroimage 61(4), 1428--1443 (2012).

    """
    # Read the surface file
    vtx, fac = read_geometry(filename_surf)
    nV = len(vtx)
    nF = len(fac)

    # compute area per face (DPF)
    facvtx = np.concatenate([vtx[fac[:, 0]], vtx[fac[:, 1]], vtx[fac[:, 2]]], axis=1)
    facvtx0 = facvtx[:, 0:6] - np.concatenate(
        [facvtx[:, 6:9], facvtx[:, 6:9]], axis=1
    )  # place 3rd vtx at origin
    cp = np.cross(facvtx0[:, 0:3], facvtx0[:, 3:6], axisa=1, axisb=1)  # cross product
    dpf = norm(cp, axis=1) / 2  # half of the norm
    print("Total area (facewise): " + str(np.sum(dpf)))

    # compute area per vertex (DPV)
    dpv = np.zeros(nV)

    # for speed, divide the dpf by 3
    dpf = dpf / 3

    # redistribute
    for f in range(nF):
        dpv[fac[f, :]] = dpv[fac[f, :]] + dpf[f]

    print("Total area (vertexwise): " + str(np.sum(dpv)))

    # save dpv
    if filename_area:
        write_morph_data(filename_area, dpv)

    # return vertex-wise surface area
    return dpv


def gradient(vtx, fac, arr_scalar, normalize=True):
    """This function computes the vertex-wise gradient of a scalar field sampled on a
    triangular mesh. The calculation is taken from [1].

    Parameters
    ----------
    vtx : ndarray
        Array of vertex coordinates.
    fac : ndarray
        Corresponding faces.
    arr_scalar : ndarray
        Scalar field values per vertex.
    normalize : bool, optional
        Normalize gradient vectors. The default is True.

    Returns
    -------
    gv : ndarray
        Vertex-wise gradient vector.
    gv_magn : ndarray
        Vertex-wise gradient magnitude.

    References
    -------
    .. [1] Mancinelli, C. et al. Gradient field estimation on triangle meshes.
    Eurographics Proceedings (2018).

    """
    # face areas and normals
    arr_a = _face_area(vtx, fac)
    arr_n = _face_normal(vtx, fac)

    # face-wise gradient
    gf_ji = arr_scalar[fac[:, 1]] - arr_scalar[fac[:, 0]]
    gf_ki = arr_scalar[fac[:, 2]] - arr_scalar[fac[:, 0]]

    v_ik = vtx[fac[:, 0], :] - vtx[fac[:, 2], :]
    v_ji = vtx[fac[:, 1], :] - vtx[fac[:, 0], :]

    # rotate
    v_ik_rot = np.cross(v_ik, arr_n)
    v_ji_rot = np.cross(v_ji, arr_n)

    gf = np.zeros_like(fac).astype(float)
    gf[:, 0] = (gf_ji * v_ik_rot[:, 0] + gf_ki * v_ji_rot[:, 0]) / (2 * arr_a)
    gf[:, 1] = (gf_ji * v_ik_rot[:, 1] + gf_ki * v_ji_rot[:, 1]) / (2 * arr_a)
    gf[:, 2] = (gf_ji * v_ik_rot[:, 2] + gf_ki * v_ji_rot[:, 2]) / (2 * arr_a)

    # vertex-wise gradient
    gv = _f2v(fac, gf, arr_a)
    gv_magn = norm(gv, axis=1)

    # normalize
    if normalize:
        gv_norm = norm(gv, axis=1)
        gv_norm[gv_norm == 0] = np.nan

        gv[:, 0] /= gv_norm
        gv[:, 1] /= gv_norm
        gv[:, 2] /= gv_norm
        pole = np.argwhere(np.isnan(gv))[:, 0]
        gv[pole, :] = 0

    return gv, gv_magn


def _face_area(v, f):
    """Helper function to compute face areas."""
    # indexed view into the vertex array
    tris = v[f]

    a = np.cross(tris[::, 1] - tris[::, 0], tris[::, 2] - tris[::, 0])
    a = norm(a, axis=1)
    a /= 2

    return a


def _face_normal(v, f):
    """Helper function to compute face-wise normals."""
    # indexed view into the vertex array
    tris = v[f]

    # calculate the normal for all triangles by taking the cross product of
    # the vectors v1-v0 and v2-v0 in each triangle
    n = np.cross(tris[::, 1] - tris[::, 0], tris[::, 2] - tris[::, 0])
    n_norm = norm(n, axis=1)

    # normalize
    n[:, 0] /= n_norm
    n[:, 1] /= n_norm
    n[:, 2] /= n_norm

    return n


def _f2v(f, gf, a):
    """Helper function to transform face- to vertex-wise expressions."""
    nv = np.max(f) + 1  # number of vertices
    nf = len(f)  # number of faces
    gv = np.zeros((nv, 3))
    magn = np.zeros(nv)
    for i in range(nf):
        gv[f[i, 0], :] += a[i] * gf[i, :]
        gv[f[i, 1], :] += a[i] * gf[i, :]
        gv[f[i, 2], :] += a[i] * gf[i, :]

        magn[f[i, 0]] += a[i]
        magn[f[i, 1]] += a[i]
        magn[f[i, 2]] += a[i]

    gv[:, 0] /= magn
    gv[:, 1] /= magn
    gv[:, 2] /= magn

    return gv


def b0_orientation(surf_in, vol_in, write_output=False, path_output="", name_output=""):
    """This function computes the angle between surface normals and B0-direction per
    vertex.

    Parameters
    ----------
    surf_in : str
        Input of surface mesh.
    vol_in : str
        Input of corresponding nifti volume.
    write_output : bool, optional
        Write out to disk (boolean). The default is False.
    path_output : str, optional
        Path where to save output. The default is "".
    name_output : str, optional
        Basename of output file. The default is "".

    Returns
    -------
    theta : ndarray
        Angle in radians.

    """
    # make subfolders
    if write_output and not os.path.exists(path_output):
        os.makedirs(path_output)

    # get hemi from surface filename
    _, hemi, _ = get_filename(surf_in)

    # load surface
    vtx, fac = read_geometry(surf_in)

    # get transformation matrix
    _, r2v = read_vox2ras_tkr(vol_in)  # ras-tkr -> voxel
    v2s = nb.load(vol_in).affine  # voxel -> scanner-ras
    m = v2s.dot(r2v)

    # apply affine transformation
    vtx = apply_affine_chunked(m, vtx)

    # get surface normals
    n = Mesh(vtx, fac).vertex_normals

    # get angle between b0 and surface normals in radians
    theta = np.arccos(np.dot(n, [0, 0, 1]))

    # write output
    if write_output:
        write_morph_data(os.path.join(path_output, hemi + "." + name_output), theta)

    return theta


def clip_surface(verts, faces, ind_keep):
    """Clip a triangle surface mesh based on a separate index list. Indices should form
    a connected region on the surface. Vertices which are not in the index list will be
    removed and the corresponding face array will be updated.

    Parameters
    ----------
    verts : np.ndarray, shape=(N,3)
        Vertex coordinates.
    faces : np.ndarray, shape=(M,3)
        Vertex indices of each triangle.
    ind_keep : list
        Index list of vertices to keep, which should be a connected region.

    Returns
    -------
    verts : np.ndarray, shape=(N,3)
        Remaining vertex coordinates.
    faces : np.ndarray, shape=(M,3)
        Vertex indices of remaining triangles.

    """
    # get new vertices
    verts = verts[ind_keep, :]

    # get new faces
    faces_keep = np.zeros(len(faces))
    faces_keep += np.in1d(faces[:, 0], ind_keep)
    faces_keep += np.in1d(faces[:, 1], ind_keep)
    faces_keep += np.in1d(faces[:, 2], ind_keep)
    faces = faces[faces_keep == 3, :]

    # reindex faces
    for i, ind in enumerate(ind_keep):
        faces[faces == ind] = i

    return verts, faces


def clip_mgh(mgh_in, label_in, mgh_out):
    """Clip mgh overlay based on label file.

    Parameters
    ----------
    mgh_in : str
        Path to input mgh file.
    label_in : str
        Path to input label file.
    mgh_out : str
        Path to output mgh file.

    Returns
    -------
    None.

    """

    arr, affine, header = read_mgh(mgh_in)
    label = read_label(label_in)
    write_mgh(mgh_out, arr[label], affine, header)
