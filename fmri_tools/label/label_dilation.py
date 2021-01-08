# -*- coding: utf-8 -*-

# external inputs
import numpy as np
from gbb.neighbor import nn_2d

# local inputs
from fmri_tools.label.label_border import label_border


def label_dilation(arr_label, adjm, n):
    """ Label dilation
    
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

    Notes
    -------
    created by Daniel Haenelt
    Date created: 20-11-2020             
    Last modified: 20-11-2020 

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
