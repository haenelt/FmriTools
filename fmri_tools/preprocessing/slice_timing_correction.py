# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb
from scipy.interpolate import InterpolatedUnivariateSpline as Interp

# local inputs
#from ..io.get_filename import get_filename
from fmri_tools.io.get_filename import get_filename


def _set_tr(img, tr):
    """Helper function to set tr in nifti header."""
    
    header = img.header.copy()
    zooms = header.get_zooms()[:3] + (tr,)
    header.set_zooms(zooms)
    
    return img.__class__(img.get_fdata().copy(), img.affine, header)


def slice_timing_correction(file_in, TR_old, TR_new, order, mb=None, 
                            manufacturer="siemens", prefix="a"):
    """Slice timing correction.

    This function performs slice timing correction of a nifti time series. For 
    interleaved slice ordering, interleaved ascending is assumed. The correction 
    is done by temporal interpolation of single voxel time series using cubic 
    interpolation. To omit extrapolation errors at the edges, the first and last 
    volumes of the time series are appended at the beginning and at the end, 
    respectively. These time points are removed again after the interpolation 
    step. The interpolated time series is sampled onto a regular grid with a 
    defined new TR. Therefore, the reference slice is always the first slice 
    acquired at t = 0. For time series acquired with multiband, the number of 
    slices has to be a multiple of the multiband factor. For interleaved slice 
    acquisition, siemens sequences start with the odd (even) slice for images 
    with odd (even) number of slices. You can switch to cmrr to consider that 
    the cmrr sequence always starts with the odd slice (first slice) in 
    interleaved mode irrespective of the number of slices.

    Parameters
    ----------
    file_in : str
        Filename of nifti time series.
    TR_old : float
        TR of time series in seconds.
    TR_new : float
        TR of slice timing corrected time series in s.
    order : str
        Slice ordering (ascending, descending, interleaved).
    mb : int
        Multiband factor.
    manufacturer : str
        Sequence type (siemens, cmrr).
    prefix : str, optional
        Prefix of output time series basename. The default is "a".

    Returns
    -------
    None.

    """

    if not mb:
        mb = 1

    if manufacturer not in ["siemens", "cmrr"]:
        raise ValueError("Unknown manufacturer!")

    # get filename
    path_file, name_file, ext_file = get_filename(file_in)

    # load data
    data = nb.load(file_in)
    nx = data.header["dim"][1]
    ny = data.header["dim"][2]
    nz = data.header["dim"][3]
    nt = data.header["dim"][4]

    # load array with appended volumes
    data_array = np.zeros((nx, ny, nz, nt + 2))
    data_array[:, :, :, 0] = data.get_fdata()[:, :, :, 0]
    data_array[:, :, :, -1] = data.get_fdata()[:, :, :, -1]
    data_array[:, :, :, 1:-1] = data.get_fdata()

    # effective number of sequentially acquired slices
    mb_package = nz/mb
    if np.mod(mb_package, 1):
        raise ValueError("Number of slices and multiband factor does not match!")
    else:
        mb_package = int(mb_package)

    # spatial order of acquired slices
    if order == "ascending":
        slice_order = np.arange(0, nz)
    elif order == "descending":
        slice_order = np.arange(nz-1, - 1, -1)
    elif order == "interleaved" and np.mod(nz, 2) or manufacturer == "cmrr":  # odd slice number
        slice_order = np.arange(0, nz, 2)
        slice_order = np.append(slice_order, np.arange(1, nz, 2))
    elif order == "interleaved" and not np.mod(nz, 2):  # even slice number
        slice_order = np.arange(1, nz, 2)
        slice_order = np.append(slice_order, np.arange(0, nz, 2))
    else:
        raise ValueError("Choose a valid slice ordering!")
    
    # temporal order of acquired slices
    if order == "interleaved":
        target = np.ceil(mb_package/2).astype(int)
        temporal_order = np.arange(0, target)
        if float(mb_package/2).is_integer():
            temporal_order2 = np.arange(0, target)
        else:
            temporal_order2 = np.arange(0, target)[:-1]            

        temporal_order = np.tile(temporal_order, mb)
        temporal_order2 = np.tile(temporal_order2, mb)
        temporal_order = np.append(temporal_order, temporal_order2 + target)
    else:
        temporal_order = np.arange(0, mb_package)
        temporal_order = np.tile(temporal_order, mb)
    
    # some prints for sanity check
    print("Spatial order of slices: "+str(slice_order))
    print("Temporal order of slices: "+str(temporal_order))
    
    # some parameters  
    TA = TR_old / mb_package  # acquisition time needed for one slice
    TT = TR_old * nt  # total acquisition time
    TR_append = np.floor(TR_old / TR_new).astype(
        int) * TR_new  # number of appended TRs in output array
    t_new = np.arange(-TR_append, TT + TR_append, TR_new)  # grid points of output array

    # temporal interpolation
    data_array_corrected = np.zeros((nx, ny, nz, len(t_new)))
    for z in range(nz):
        print("Slice timing correction for slice: " + str(z + 1) + "/" + str(nz))
        for x in range(nx):
            for y in range(ny):
                t = np.arange(temporal_order[z] * TA - TR_old, temporal_order[z] * TA + (nt + 1) * TR_old, TR_old)
                cubic_interper = Interp(t, data_array[x, y, slice_order[z], :], k=3)
                data_array_corrected[x, y, slice_order[z], :] = cubic_interper(t_new)

    # delete appended volumes
    vols_keep1 = t_new >= 0
    vols_keep2 = t_new < TT
    vols_keep = vols_keep1 * vols_keep2
    data_array_corrected = data_array_corrected[:, :, :, vols_keep]

    # clean corrected array
    data_min = np.min(data.get_fdata())
    data_max = np.max(data.get_fdata())
    data_array_corrected[np.isnan(data_array_corrected)] = 0
    data_array_corrected[data_array_corrected < data_min] = data_min
    data_array_corrected[data_array_corrected > data_max] = data_max

    # update data header
    data.header["dim"][4] = np.shape(data_array_corrected)[3]
    data.header["datatype"] = 16

    # write output
    output = nb.Nifti1Image(data_array_corrected, data.affine, data.header)
    output = _set_tr(output, TR_new)
    nb.save(output, os.path.join(path_file, prefix + name_file + ext_file))
