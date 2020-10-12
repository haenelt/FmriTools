# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys

# external inputs
import numpy as np
import nibabel as nb
from scipy.interpolate import InterpolatedUnivariateSpline as Interp

# local inputs
from fmri_tools.io.get_filename import get_filename


def slice_timing_correction(input, TR_old, TR_new, order, prefix="a"):    
    """
    This function performs slice timing correction of a nifti time series. For 
    interleaved slice ordering, interleaved ascending is assumed. The correction 
    is done by temporal interpolation of single voxel time series using cubic 
    interpolation. To omit extrapolation errors at the edges, the first and last 
    volumes of the time series are appended at the beginning and at the end, 
    respectively. These time points are removed again after the interpolation 
    step. The interpolated time series is sampled onto a regular grid with a 
    defined new TR. Therefore, the reference slice is always the first slice 
    acquired at t = 0. Only for writing the new TR in the header of the output 
    time series, AFNI has to be included in the search path.
    Inputs:
        *input: filename of nifti time series.
        *TR_old: TR of time series in seconds.
        *TR_new: TR of slice timing corrected time series in s.
        *order: slice ordering (ascending, descending, interleaved).
        *prefix: prefix of output time series basename.
            
    created by Daniel Haenelt
    Date created: 11-03-2019
    Last modified: 12-10-2020
    """

    # get filename
    path_file, name_file, ext_file = get_filename(input)

    # load data
    data = nb.load(input)
    nx = data.header["dim"][1]
    ny = data.header["dim"][2]
    nz = data.header["dim"][3]
    nt = data.header["dim"][4]

    # load array with appended volumes
    data_array = np.zeros((nx, ny, nz, nt+2))
    data_array[:,:,:,0] = data.get_fdata()[:,:,:,0]
    data_array[:,:,:,-1] = data.get_fdata()[:,:,:,-1]
    data_array[:,:,:,1:-1] = data.get_fdata()

    # get slice order
    if order == "ascending":
        slice_order = np.arange(0, nz)
    elif order == "descending":
        slice_order = np.arange(nz-1,-1,-1)
    elif order == "interleaved" and np.mod(nz,2): # odd slice number
        slice_order = np.arange(0,nz,2)
        slice_order = np.append(slice_order, np.arange(1,nz,2))
    elif order == "interleaved" and not np.mod(nz,2): # even slice number
        slice_order = np.arange(1,nz,2)
        slice_order = np.append(slice_order, np.arange(0,nz,2))    
    else:
        sys.exit("Choose a valid slice ordering!")

    # some parameters  
    TA = TR_old / nz # acquisition time needed for one slice
    TT = TR_old * nt # total acquisition time
    TR_append = np.floor(TR_old/TR_new).astype(int) * TR_new # number of appended TRs in output array
    t_new = np.arange(-TR_append, TT+TR_append, TR_new) # grid points of output array

    # temporal interpolation
    data_array_corrected = np.zeros((nx,ny,nz,len(t_new)))
    for z in range(nz):
        print("Slice timing correction for slice: "+str(z+1)+"/"+str(nz))
        for x in range(nx):
            for y in range(ny):
                t = np.arange(z*TA-TR_old, z*TA+(nt+1)*TR_old, TR_old)
                cubic_interper = Interp(t, data_array[x,y,slice_order[z],:], k=3)
                data_array_corrected[x,y,slice_order[z],:] = cubic_interper(t_new)

    # delete appended volumes
    vols_keep1 = t_new >= 0
    vols_keep2 = t_new < TT
    vols_keep = vols_keep1 * vols_keep2
    data_array_corrected = data_array_corrected[:,:,:,vols_keep]

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
    nb.save(output, os.path.join(path_file,prefix+name_file+ext_file))
    
    # change TR in header
    os.system("3drefit " + \
              "-TR " + str(TR_new) + " " + \
              os.path.join(path_file,prefix+name_file+ext_file))
