# -*- coding: utf-8 -*-

# external inputs
import numpy as np
from gbb.neighbor import nn_2d

# internal inputs
from ..io.affine import vox2ras_tkr
from ..utils.apply_affine_chunked import apply_affine_chunked

__all__ = ['label_border', 'label_dilation', 'label_erosion', 'roi_fov', 'roi_sphere']


def label_border(arr_label, adjm):
    """Label border.

    This function returns border vertex indices from an input array containing
    vertex indices of a freesurfer label.

    Parameters
    ----------
    arr_label : ndarray
        1D array of label indices.
    adjm : ndarray
        Adjacency matrix.

    Returns
    -------
    border : ndarray
        1D array of border indices.

    """

    # label array as set
    arr_label_set = set(arr_label)

    border = []
    for i in arr_label:

        # get nearest neighbors
        nn = nn_2d(i, adjm, 0)

        # check if all neighbors are within the label
        if not set(nn).issubset(arr_label_set):
            border.append(i)

    return np.array(border)


def label_dilation(arr_label, adjm, n):
    """Label dilation.

    This function dilates a labeled region of interest which is defined as a 1D
    array of triangular mesh indices. Dilation is done by adding the one-ring
    neighborhood of all border vertices to the label array. This can be done
    iteratively.

    Parameters
    ----------
    arr_label : ndarray
        1D array of label indices.
    adjm : ndarray
        Adjacency matrix.
    n : int
        Number of dilation iterations.

    Returns
    -------
    arr_label : ndarray
        1D array of dilated label indices.

    """

    arr_dilate = []
    for _ in range(n):
        border = label_border(arr_label, adjm)
        for j in border:
            nn = nn_2d(j, adjm, 0)
            arr_dilate.extend(nn)
        arr_label = np.append(arr_label, arr_dilate)
        arr_label = np.unique(arr_label)
        arr_label = np.sort(arr_label)
    return arr_label


def label_erosion(arr_label, adjm, n):
    """Label erosion.

    This function erodes a labeled region of interest which is defined as a 1D
    array of triangular mesh indices. Erosion is done by removing all border
    indices from the label array. This can be done iteratively.

    Parameters
    ----------
    arr_label : ndarray
        1D array of label indices.
    adjm : ndarray
        Adjacency matrix.
    n : int
        Number of erosion iterations.

    Returns
    -------
    arr_label : ndarray
        1D array of eroded label indices.

    """

    for _ in range(n):
        border = label_border(arr_label, adjm)
        tmp = np.in1d(arr_label, border)
        arr_label = arr_label[tmp != 1]
    return arr_label


def roi_fov(vtx, vol_dims, vol_ds):
    """ROI from image FOV.

    This function creates a ROI label for all vertex indices within an imaging FOV. Nans
    in the vertex array are excluded.

    Parameters
    ----------
    vtx : ndarray
        Vertex array
    vol_dims : tuple
        Tuple containing volume dimensions in x-, y- and z-direction.
    vol_ds : tuple
        Tuple containing voxel sizes in x-, y- and z-direction.

    Returns
    -------
    arr_label : ndarray
        1D array of roi indices.

    """

    _, r2v = vox2ras_tkr(vol_dims, vol_ds)  # affine transformation to voxel space
    vtx_voxel = apply_affine_chunked(r2v, vtx)  # apply transformation to vertex array

    # mask vertices within volume dimensions
    arr_label = np.arange(len(vtx), dtype=int)
    for i, v in enumerate(vol_dims):
        arr_label[np.isnan(vtx_voxel[:, 0])] = -1
        arr_label[np.isnan(vtx_voxel[:, 1])] = -1
        arr_label[np.isnan(vtx_voxel[:, 2])] = -1
        arr_label[vtx_voxel[:, i] < 0] = -1
        arr_label[vtx_voxel[:, i] > v - 1] = -1

    return arr_label[arr_label != -1]


def roi_sphere(vtx, ind, radius):
    """Spherical ROI.

    This function creates a ROI label for all vertex indices within a 3D sphere. Nans in
    the vertex array are excluded automatically.

    Parameters
    ----------
    vtx : ndarray
        Vertex array
    ind : int
        Vertex index of the center of the sphere.
    radius : float
        Radius of the sphere.

    Returns
    -------
    arr_label : ndarray
        1D array of roi indices.

    """

    x_diff = vtx[:, 0] - vtx[ind, 0]
    y_diff = vtx[:, 1] - vtx[ind, 1]
    z_diff = vtx[:, 2] - vtx[ind, 2]
    distance = np.sqrt(x_diff**2 + y_diff**2 + z_diff**2)

    arr_label = np.arange(len(vtx))
    return arr_label[distance <= radius]
