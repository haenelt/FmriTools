# -*- coding: utf-8 -*-

# external inputs
import numpy as np
from gbb.neighbor import nn_2d


def label_border(arr_label, adjm):
    """ Label border
    
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

    Notes
    -------
    created by Daniel Haenelt
    Date created: 19-11-2020             
    Last modified: 19-11-2020 

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
