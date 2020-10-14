# -*- coding: utf-8 -*-

# python standard library inputs
import math

# external inputs
import numpy as np


def mask_pattern_2d(nx, ny, a, b, alpha):
    """ Mask pattern 2D

    This function computes a mask to occlude parts of the columnar pattern. The 
    mask has elliptical shape and can be adjusted by its major and minor axes 
    and a rotational angle.    

    Parameters
    ----------
    nx : int
        Array size in x-direction.
    ny : int
        Array size in y-direction.
    a : float
        Major axis of ellipse.
    b : float
        Minor axis of ellipse.
    alpha : float
        Rotation in degrees.

    Returns
    -------
    mask : ndarray
        Binary array containing mask.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 12-04-2019
    Last modified: 12-10-2020
    
    """
    
    # compute meshgrid
    x = np.linspace(-nx/2,nx/2,nx)
    y = np.linspace(-ny/2,ny/2,ny)
    x_mesh, y_mesh = np.meshgrid(x,y)

    # rotate meshgrid
    x_rot = x_mesh*np.cos(math.radians(alpha)) - y_mesh*np.sin(math.radians(alpha))
    y_rot = x_mesh*np.sin(math.radians(alpha)) + y_mesh*np.cos(math.radians(alpha))

    # compute ellipse
    mask = x_rot**2 / a**2 + y_rot**2 / b**2

    # mask ellipse
    mask[mask > 1] = 0
    mask[mask != 0] = 1

    return mask


def mask_pattern_1d(n, a, b):
    """ Mask pattern 1D
    
    This function computes a mask to occlude parts of the columnar pattern.

    Parameters
    ----------
    n : int
        Array size.
    a : float
        Mask left side.
    b : float
        Mask right side.

    Returns
    -------
    mask : ndarray
        Binary array containing mask.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 12-04-2019
    Last modified: 12-10-2020

    """
    
    # compute meshgrid
    mask = np.linspace(-n/2,n/2,n)

    # mask ellipse
    mask[mask < -a] = 0
    mask[mask > b] = 0
    mask[mask != 0] = 1

    return mask
