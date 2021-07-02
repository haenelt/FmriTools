# -*- coding: utf-8 -*-

# external inputs
import numpy as np
from numpy.linalg import norm


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
    """ Helper function to transform face- to vertex-wise expressions."""

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


def gradient(vtx, fac, arr_scalar, normalize=True):
    """Gradient.
    
    This function computes the vertex-wise gradient of a scalar field sampled
    on a triangular mesh. The calculation is taken from [1].

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
