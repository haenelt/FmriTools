# -*- coding: utf-8 -*-

# external inputs
import numpy as np
from gbb.neighbor import nn_2d

__all__ = ['label_border', 'label_dilation', 'label_erosion']


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
    for i in range(n):

        # get label border
        border = label_border(arr_label, adjm)

        # dilate
        for j in border:
            nn = nn_2d(j, adjm, 0)
            arr_dilate.extend(nn)

        # update label indices
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

    for i in range(n):
        # get label border
        border = label_border(arr_label, adjm)

        # update label indices
        tmp = np.in1d(arr_label, border)
        arr_label = arr_label[tmp == False]

    return arr_label
