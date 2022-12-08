# -*- coding: utf-8 -*-

# external inputs
import numpy as np
import nibabel as nb

# local inputs
from ..io.affine import read_vox2ras_tkr
from ..utils.apply_affine_chunked import apply_affine_chunked
from ..utils.interpolation import nn_interpolation3d, linear_interpolation3d

# linear and nearest neighbor sampling functions
_sampler = {"linear": linear_interpolation3d, "nearest": nn_interpolation3d}


def map_timeseries(vtx, file_timeseries, interpolation="linear"):
    """Map time series data onto a surface mesh. A 2D array is returned which contains
    vertex-wise sampled data for each time point in separate columns. All vertices
    outside the time series volume are set to nan.

    Parameters
    ----------
    vtx : np.ndarray, shape=(N,3)
        Vertex-wise array.
    file_timeseries : str
        File name of fMRI time series.
    interpolation : str, optional (linear | nearest)
        Interpolation method (linear or nearest neighbor interpolation).

    Returns
    -------
    arr_sampled : np.ndarray, shape=(N,T)
        Vertex-wise sampled time series.

    """

    arr = nb.load(file_timeseries).get_fdata()
    nx, ny, nz, nt = np.shape(arr)
    _, ras2vox = read_vox2ras_tkr(file_timeseries)
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
            vtx_vox[:, 0], vtx_vox[:, 1], vtx_vox[:, 2], arr[:, :, :, i]
        )

    return arr_sampled
