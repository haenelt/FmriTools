# -*- coding: utf-8 -*-

# external inputs
import numpy as np


def apply_affine_chunked(aff, pts, chunk_size=10000):
    """ Apply an affine matrix to points in chunks.
    
    This function is a copy of the routine `apply_affine` from the `affines`
    module of the `nibabel` package. The only difference is that this function
    applies the affine transformation in chunks to prevent memory errors when 
    working with large arrays. More information about this function can be found 
    in the docstring of the aforementioned nibabel function.

    Parameters
    ----------
    aff : (N, N) np.ndarray
        Homogenous affine, for 3D points, will be 4 by 4. Contrary to first
        appearance, the affine will be applied on the left of `pts`.
    pts : (..., N-1) np.ndarray
        Points, where the last dimension contains the coordinates of each
        point. For 3D, the last dimension will be length 3.

    Returns
    -------
    (..., N-1) np.ndarray
        Transformed points.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 11-03-2021
    Last modified: 11-03-2021

    """
    
    aff = np.asarray(aff)
    pts = np.asarray(pts)
    shape = pts.shape
    pts = pts.reshape((-1, shape[-1]))
    # rzs == rotations, zooms, shears
    rzs = aff[:-1, :-1]
    trans = aff[:-1, -1]

    # chunk intervals
    chunk = np.arange(0, len(pts), chunk_size)
    chunk[-1] = len(pts)
    if chunk[0] == 0:
        chunk = np.delete(chunk, 0)

    # apply affine transformation piecewise
    j1 = 0
    res = np.zeros_like(pts)
    for _, j2 in enumerate(chunk):
        res[j1:j2,:] = np.dot(pts[j1:j2,:], rzs.T) + trans[None, :]
        j1 = j2

    return res.reshape(shape)
