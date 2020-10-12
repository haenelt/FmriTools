# -*- coding: utf-8 -*-

# python standard library inputs
import math

# external inputs
import numpy as np


def mask_pattern_2d(nx, ny, a, b, alpha):
    """
    This function computes a mask to occlude parts of the columnar pattern. The 
    mask has elliptical shape and can be adjusted by its major and minor axes 
    and a rotational angle.
    Inputs:
        *nx: array size in x-direction.
        *ny: array size in y-direction.
        *a: major axis of ellipse.
        *b: minor axis of ellipse.
        *alpha: rotation in degrees.
    Outputs:
        *mask: binary array containing mask.
        
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
    """
    This function computes a mask to occlude parts of the columnar pattern.
    Inputs:
        *nx: array size.
        *a: mask left side.
        *b: mask right side.
    Outputs:
        *ellipse: binary array containing mask.
        
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
