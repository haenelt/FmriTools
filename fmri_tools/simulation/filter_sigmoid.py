# -*- coding: utf-8 -*-

# external inputs
import numpy as np
    

def filter_sigmoid(x, alpha=4):
    """ Filter sigmoid

    Nonlinear sigmoid filter to adjust sharpness of ocular dominance. If 
    sharpness factor is set to zero, nothing is changed on the input image. The 
    sigmoidal filter outputs an array with intensity range [-1,1].    

    Parameters
    ----------
    x : ndarray
        Input image.
    alpha : float, optional
        Sharpness factor. The default is 4.

    Returns
    -------
    ndarray
        Filtered image.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 29-11-2018             
    Last modified: 12-10-2020
    
    """
    
    if alpha:            
        #return 1 / (1 + np.exp(-alpha*x))
        return 2 / ( 1 + np.exp(-alpha*x) ) -1
    else:
        return x
