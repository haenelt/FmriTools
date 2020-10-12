# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb

# local inputs
from fmri_tools.processing.demean_time_series import demean_time_series


"""
Average cycle periods

The purpose of the following script is to average a time points within a time 
series depending on the stimulus cycle for display purposes.

created by Daniel Haenelt
Date created: 22-05-2020
Last modified: 12-10-2020
"""

# input
file_in = "/data/pt_01880/Experiment4_PSF/p6/psf/GE_EPI1/multipol_2/uadata.nii"

# parameters
path_output = "/data/pt_01880/Experiment4_PSF/p6/psf/GE_EPI1/multipol_2"
name_output = "test" # basename of output file
TR = 3 # repetition time in s
dummy = 12 # number of skipped volumes at the beginning and at the end
period = 48 # cycle period in s
ncycle = 10.5 # number of cycles (fractional cycles will be skipped)
demean = True # detrend data

# do not edit below

# skipped time points in vols
skip = int(dummy/TR)
skip_cycle = int(np.mod(ncycle,1)*period/TR)
full_cycle = int(np.floor(ncycle))

# load data
data = nb.load(file_in)
data_array = data.get_fdata()

# get array dimensions
xdim = data.header["dim"][1]
ydim = data.header["dim"][2]
zdim = data.header["dim"][3]

# exclude dummy volumes
data_array = data_array[:,:,:,skip+skip_cycle:-skip]

# merge periods together
merge_array = np.zeros((xdim,ydim,zdim,full_cycle))

c = 0
for i in range(np.shape(data_array)[3]):
    merge_array[:,:,:,c] += data_array[:,:,:,i]
    c += 1
    if c == full_cycle:
        c = 0

merge_array /= full_cycle

# write output
data.header["dim"][4] = int(period/TR)
output = nb.Nifti1Image(merge_array, data.affine, data.header)
nb.save(output,os.path.join(path_output,name_output+".nii"))

if demean:
    demean_time_series(output, path_output, name_output, write_output=True)
