# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb
from sh import gunzip
from scipy.interpolate import InterpolatedUnivariateSpline as Interp

# local inputs
from fmri_tools.io.get_filename import get_filename

    
def regrid_time_series(input, path_output, TR_old, TR_new, t_start=0):
    """
    This function interpolates the time series onto a new time grid using cubic 
    interpolation. Only for writing the new TR in the header of the output time 
    series, AFNI has to be included in the search path.
    Inputs:
        *input: time series filename.
        *path_output: path where output is written.
        *TR_old: TR of time series in s.
        *TR_new: TR of regridded time series in s.
        *t_start: shift time series in s (t_start >= 0 and <= TR_old).
        
    created by Daniel Haenelt
    Date created: 19-02-2020           
    Last modified: 12-10-2020
    """

    # get filename
    _, name_input, ext_input = get_filename(input)    
    
    # print to console
    print("time series regridding for: "+name_input)
    
    # load data
    data = nb.load(input)
    nx = data.header["dim"][1]
    ny = data.header["dim"][2]
    nz = data.header["dim"][3]
    nt = data.header["dim"][4]
    
    # get time grid
    TT = TR_old * nt # total acquisition time
    TR_append = np.floor(t_start/TR_old + 1).astype(int) * TR_old # number of appended TRs in input array
    
    # input grid
    t_old = np.arange(-TR_append, TT + TR_append, TR_old) + t_start
    
    # output grid
    t_new = np.arange(0, TT + TR_append, TR_new)
    t_new_append = np.flip(np.arange(0,-TR_append,-TR_new)[1:])
    if not len(t_new_append):
        t_new_append = -TR_new    
    t_new = np.append(t_new_append, t_new)
    
    # load array with appended volumes
    n_append = int(TR_append/TR_old)
    data_array = np.zeros((nx, ny, nz, nt+2*n_append))
    data_array[:,:,:,n_append:-n_append] = data.get_fdata()
    for i in range(n_append):
        data_array[:,:,:,i] = data.get_fdata()[:,:,:,0]
        data_array[:,:,:,-(i+1)] = data.get_fdata()[:,:,:,-1]
    
    # temporal interpolation
    data_array_regrid = np.zeros((nx,ny,nz,len(t_new)))
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                cubic_interper = Interp(t_old, data_array[x,y,z,:], k=3)
                data_array_regrid[x,y,z,:] = cubic_interper(t_new)
    
    # delete appended volumes
    vols_keep1 = t_new >= 0
    vols_keep2 = t_new < TT
    vols_keep = vols_keep1 * vols_keep2
    data_array_regrid = data_array_regrid[:,:,:,vols_keep]
    
    # clean corrected array
    data_array_regrid[np.isnan(data_array_regrid)] = 0
    data_array_regrid[data_array_regrid < 0] = 0
    
    # update data header
    data.header["dim"][4] = np.shape(data_array_regrid)[3]
    data.header["datatype"] = 16
    
    # write output
    output = nb.Nifti1Image(data_array_regrid, data.affine, data.header)
    nb.save(output, os.path.join(path_output,name_input+"_upsampled"+ext_input))
    
    # change TR in header
    os.system("3drefit " + \
              "-TR " + str(TR_new) + " " + \
              os.path.join(path_output,name_input+"_upsampled"+ext_input))


def regrid_time_series_afni(input, n=2):
    """
    This function upsamples a time series using the afni function 3dUpsample. 
    Before running the function, set the afni environment by calling AFNI in the 
    terminal. Output is an upsampled nifti time series.
    Inputs:
        *input: time series filename.
        *n: upsampling factor.
        
    created by Daniel Haenelt
    Date created: 20-09-2019           
    Last modified: 12-10-2020
    """
    
    clean_unzip = 0
    if os.path.splitext(input)[1] == ".gz":
        gunzip(input)
        clean_unzip = 1
        input = os.path.splitext(input)[0]
        
    # prepare path and filename
    path_file = os.path.dirname(input)
    name_file = os.path.splitext(os.path.basename(input))[0]

    # upsample vaso and bold time series
    os.system("3dUpsample -overwrite -datum short " + \
              "-prefix " + os.path.join(path_file,name_file + "_upsampled.nii") + \
              " -n " + str(n) + " -input " + input)
    
    if clean_unzip:
        os.remove(input)
        