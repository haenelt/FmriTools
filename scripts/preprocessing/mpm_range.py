# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb

# local inputs
from fmri_tools.io.get_filename import get_filename


"""
MPM scaling

This scripts restricts values in a nifti volume to a defined data range. Ceiling 
outside are set to the set maximum or minimum value. Optionally, the range is 
inverted.

created by Daniel Haenelt
Date created: 03-03-2020
Last modified: 12-10-2020
"""

# input
input = [
    "/data/pt_01880/Experiment3_Stripes/p1/mpm/pd_kp_mtflash3d_v1ax_0p5_0008/Results/s2899875-124739-00001-00352-1_PD_gnlcorr.nii",
    ]

val_min = 0 # if not None
val_max = 100 # if not None
val_invert = False # if True

# do not edit below

for i in range(len(input)):
    
    # get filename
    path_file, name_file, ext_file = get_filename(input[i])
    
    # load data
    data = nb.load(input[i])
    data_array = data.get_fdata()
    
    # scale
    if val_min is not None:
        data_array[data_array < val_min] = val_min
        
    if val_max is not None:
        data_array[data_array > val_max] = val_max

    # invert    
    if val_invert:
        data_array = np.max(data_array) - data_array
    
    # write output
    output = nb.Nifti1Image(data_array, data.affine, data.header)
    if val_invert:
        file_out =  os.path.join(path_file, name_file + "_scaled_inverted" + ext_file)
    else:
        file_out =  os.path.join(path_file, name_file + "_scaled" + ext_file)
    
    nb.save(output, file_out)
