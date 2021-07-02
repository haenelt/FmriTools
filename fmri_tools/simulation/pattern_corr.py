# -*- coding: utf-8 -*-

# external inputs
import numpy as np
from scipy.stats import pearsonr


def pattern_corr(data_array1, data_array2):
    """Pattern corr.

    This function calculates the pearson correlation between two arrays of the 
    same size.    

    Parameters
    ----------
    data_array1 : ndarray
        Input array1.
    data_array2 : ndarray
        Input array2.

    Returns
    -------
    tuple
        Pearson coefficient
    
    """
    
    # reshape data arrays into a vector
    x = np.reshape(data_array1, (np.size(data_array1), 1))
    y = np.reshape(data_array2, (np.size(data_array2), 1))
    
    return pearsonr(x, y)
