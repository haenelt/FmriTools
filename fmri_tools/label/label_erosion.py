# -*- coding: utf-8 -*-

# external inputs
import numpy as np

# local inputs
from fmri_tools.label.label_border import label_border


def label_erosion(arr_label, adjm, n):
    """ Label erosion
    
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
    
    Notes
    -------
    created by Daniel Haenelt
    Date created: 20-11-2020             
    Last modified: 20-11-2020 

    """
    
    for i in range(n):
        
        # get label border
        border = label_border(arr_label, adjm)
        
        # update label indices
        tmp = np.in1d(arr_label, border)
        arr_label = arr_label[tmp == False]
    
    return arr_label
