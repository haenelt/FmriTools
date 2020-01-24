"""
Venous mask from phase data

A binary mask of venous voxels is created from time series phase data. First the data is scaled to
[-pi; +pi]. Then, the time series standard deviation is taken. Venous voxels are classified by 
applying a threshold to identify voxels with large dispersion over time (phase_threshold).

created by Daniel Haenelt
Date created: 24-01-2020             
Last modified: 24-01-2020
"""
import os
import datetime
import numpy as np
import nibabel as nb

# input data
phase_input = "/data/pt_01880/test_data/data/udata_phase_unwrap.nii"

# parameters
phase_threshold = 0.05 # in rad

""" do not edit below """

# prepare path and filename
path = os.path.split(phase_input)[0]
file = os.path.split(phase_input)[1]

# output folder is taken from the first entry of the input list and set into 
path_output = os.path.join(path,"vein","native")
if not os.path.exists(path_output):
    os.makedirs(path_output)

# get image header information
phase_img = nb.load(phase_input)
phase_img.header["dim"][0] = 3
phase_img.header["dim"][4] = 1
header = phase_img.header
affine = phase_img.affine

phase_array = phase_img.get_fdata()

# normalize in range [-pi, pi]
phase_array = 2 * ( phase_array - np.min(phase_array) ) / ( np.max(phase_array) - np.min(phase_array) ) - 1
phase_array *= np.pi

# get standard deviation
phase_array = np.std(phase_array, axis=3)

# get vein mask
mask_array = phase_array.copy()
mask_array[mask_array > phase_threshold] = 0
mask_array[mask_array != 0] = 1
mask_array -= 1
mask_array = np.abs(mask_array)

output = nb.Nifti1Image(phase_array, affine, header)
fileOUT = os.path.join(path_output,"phase_dispersion.nii")
nb.save(output,fileOUT)

output = nb.Nifti1Image(mask_array, affine, header)
fileOUT = os.path.join(path_output,"vein.nii")
nb.save(output,fileOUT)

# write log
fileID = open(os.path.join(path_output,"vein_info.txt"),"a")
fileID.write("script executed: "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
fileID.write("phase_threshold: "+str(phase_threshold)+"\n")
fileID.close()