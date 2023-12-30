# -*- coding: utf-8 -*-

import numpy as np

from ..io.affine import vox2ras_tkr
from ..registration.transform import apply_affine_chunked
from ..utils.interpolation import linear_interpolation3d, nn_interpolation3d

# linear and nearest neighbor sampling functions
_sampler = {"linear": linear_interpolation3d, "nearest": nn_interpolation3d}


def map_timeseries(vtx, arr_timeseries, dims, ds, interpolation="linear"):
    """Map time series data onto a surface mesh. A 2D array is returned which contains
    vertex-wise sampled data for each time point in separate columns. All vertices
    outside the time series volume are set to nan.

    Parameters
    ----------
    vtx : np.ndarray, shape=(N,3)
        Vertex-wise array.
    arr_timeseries : np.ndarray, shape=(X,Y,Z,T)
        4D array of fMRI time series.
    dims : tuple
        Tuple containing volume dimensions in x-, y- and z-direction.
    ds : tuple
        Tuple containing voxel sizes in x-, y- and z-direction.
    interpolation : str, optional (linear | nearest)
        Interpolation method (linear or nearest neighbor interpolation).

    Returns
    -------
    arr_sampled : np.ndarray, shape=(N,T)
        Vertex-wise sampled time series.

    """

    nx, ny, nz, nt = np.shape(arr_timeseries)
    _, ras2vox = vox2ras_tkr(dims, ds)
    vtx_vox = apply_affine_chunked(ras2vox, vtx)

    # exclude nans and vertices outside of the volume
    mask = np.ones(len(vtx), dtype=bool)
    for i, n in enumerate((nx, ny, nz)):
        mask[np.isnan(vtx_vox[:, i])] = 0
        mask[vtx_vox[:, i] < 0] = 0
        mask[vtx_vox[:, i] > n - 1] = 0
    vtx_vox = vtx_vox[mask == 1, :]

    # initialize resulting array
    arr_sampled = np.empty((len(vtx), nt))
    arr_sampled[:] = np.nan

    for i in range(nt):
        arr_sampled[mask == 1, i] = _sampler[interpolation](
            vtx_vox[:, 0], vtx_vox[:, 1], vtx_vox[:, 2], arr_timeseries[:, :, :, i]
        )

    return arr_sampled
