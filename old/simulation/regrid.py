# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import griddata


def regrid_2d(data_array, nx, ny):
    """Regrid 2D.

    This function loads a two-dimensional numpy array and interpolates it to a
    new grid array with defined sizes using nearest neighbour interpolation.

    Parameters
    ----------
    data_array : ndarray
        Input array.
    nx : int
        Resolution of output array in x-direction.
    ny : int
        Resolution of output array in y-direction.

    Returns
    -------
    data_array_new : ndarray
        Output array.

    """

    # size of input array
    data_size = np.shape(data_array)

    # grid of input array
    x = np.arange(0, data_size[0])
    y = np.arange(0, data_size[1])
    xgrid, ygrid = np.meshgrid(x, y)
    xgrid = np.reshape(xgrid, np.size(xgrid))
    ygrid = np.reshape(ygrid, np.size(ygrid))
    xi_old = np.stack((xgrid, ygrid), 1)

    # grid of output array
    x_new = np.linspace(0, data_size[0], nx)
    y_new = np.linspace(0, data_size[1], ny)
    xgrid_new, ygrid_new = np.meshgrid(x_new, y_new)
    xgrid_new = np.reshape(xgrid_new, np.size(xgrid_new))
    ygrid_new = np.reshape(ygrid_new, np.size(ygrid_new))
    # xi_new = np.floor(np.stack((xgrid_new,ygrid_new),1))
    xi_new = np.stack((xgrid_new, ygrid_new), 1)

    # values of input array
    data_array = np.reshape(data_array, np.size(data_array))

    # grid values from old grid to new grid
    data_array_new = griddata(xi_old, data_array, xi_new, method="nearest")
    data_array_new = np.reshape(data_array_new, (nx, ny))

    return data_array_new


def regrid_1d(data_array, n):
    """Regrid 1D.

    This function loads a one-dimensional numpy array and interpolates it to a
    new grid array with defined sizes using nearest neighbour interpolation.

    Parameters
    ----------
    data_array : ndarray
        Input array.
    n : int
        Resolution of output array.

    Returns
    -------
    data_array_new : ndarray
        Array values interpolated to new grid.

    """

    # size of input array
    data_size = np.shape(data_array)

    # grid of input array
    xi_old = np.arange(0, data_size[0])

    # grid of output array
    # xi_new = np.floor(np.linspace(0,data_size[0],N))
    xi_new = np.linspace(0, data_size[0], n)

    # grid values from old grid to new grid
    data_array_new = griddata(xi_old, data_array, xi_new, method="nearest")

    return data_array_new
