# -*- coding: utf-8 -*-

# external inputs
import numpy as np

__all__ = ["linear_interpolation3d", "nn_interpolation3d"]


def linear_interpolation3d(x, y, z, arr_c):
    """Apply a linear interpolation of values in a 3D volume to an array of coordinates.

    Parameters
    ----------
    x : (N,) np.ndarray
        x-coordinates in voxel space.
    y : (N,) np.ndarray
        y-coordinates in voxel space.
    z : (N,) np.ndarray
        z-coordinates in voxel space.
    arr_c : (U,V,W) np.ndarray
        3D array with input values.
    Returns
    -------
    c : (N,) np.ndarray
        Interpolated values for [x,y,z].

    """

    # corner points
    x0 = np.floor(x).astype(int)
    x1 = np.ceil(x).astype(int)
    y0 = np.floor(y).astype(int)
    y1 = np.ceil(y).astype(int)
    z0 = np.floor(z).astype(int)
    z1 = np.ceil(z).astype(int)

    # distances to corner points
    xd = [_careful_divide(x[i], x0[i], x1[i]) for i, _ in enumerate(x)]
    yd = [_careful_divide(y[i], y0[i], y1[i]) for i, _ in enumerate(y)]
    zd = [_careful_divide(z[i], z0[i], z1[i]) for i, _ in enumerate(z)]

    xd = np.asarray(xd)
    yd = np.asarray(yd)
    zd = np.asarray(zd)

    # corner values
    c000 = arr_c[x0, y0, z0]
    c001 = arr_c[x0, y0, z1]
    c010 = arr_c[x0, y1, z0]
    c011 = arr_c[x0, y1, z1]
    c100 = arr_c[x1, y0, z0]
    c101 = arr_c[x1, y0, z1]
    c110 = arr_c[x1, y1, z0]
    c111 = arr_c[x1, y1, z1]

    # interpolation along x-axis
    c00 = c000 * (1 - xd) + c100 * xd
    c01 = c001 * (1 - xd) + c101 * xd
    c10 = c010 * (1 - xd) + c110 * xd
    c11 = c011 * (1 - xd) + c111 * xd

    # interpolation along y-axis
    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd

    # interpolation along z-axis
    c = c0 * (1 - zd) + c1 * zd

    return c


def nn_interpolation3d(x, y, z, arr_c):
    """Apply a nearest neighbor interpolation of values in a 3D volume to an array of
    coordinates.

    Parameters
    ----------
    x : (N,) np.ndarray
        x-coordinates in voxel space.
    y : (N,) np.ndarray
        y-coordinates in voxel space.
    z : (N,) np.ndarray
        z-coordinates in voxel space.
    arr_c : (U,V,W) np.ndarray
        3D array with input values.
    Returns
    -------
    c : (N,) np.ndarray
        Interpolated values for [x,y,z].

    """

    # get nearest neighbour grid points
    x0 = np.round(x).astype(int)
    y0 = np.round(y).astype(int)
    z0 = np.round(z).astype(int)

    c = arr_c[x0, y0, z0]

    return c


def _careful_divide(v, v0, v1):
    """Only divide if v0 and v1 are different from each other."""

    return (v - v0) / (v1 - v0) if v1 != v0 else v
