"""
Venous mask from phase data

A binary mask of venous voxels is created from a set of time series. First, a baseline correction of
each time series is applied if not done before (i.e., if no file with prefix b is found). Venous 
voxels are classified in the following way: Mean epi and mean tSNR across all time series are 
thresholded using arbitrary thresholds (<epi_threshold>, <tsnr_threshold>) considering only low mean
and low tsnr voxels. The thresholds have to be adjusted individually. The final mask is created by 
intersecting both masks. Both mean epi and mean tsnr are saved for single runs and the whole 
session. Additionally, a text file is written containing the used parameters. The method is taken 
from Kashyap et al., 2017.

created by Daniel Haenelt
Date created: 24-01-2020             
Last modified: 24-01-2020
"""
import os
import datetime
import numpy as np
import nibabel as nb

# input data
phase_input = "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/Run_1/udata.nii",

# parameters
phase_threshold = 0.2

""" do not edit below """

# prepare path and filename
path = []
file = []
for i in range(len(img_input)):
    path.append(os.path.split(img_input[i])[0])
    file.append(os.path.split(img_input[i])[1])

# output folder is taken from the first entry of the input list and set into 
if len(img_input) > 1:
    path_output = os.path.join(os.path.dirname(path[0]),"vein","native")
else:
    path_output = os.path.join(path[0],"vein","native")
if not os.path.exists(path_output):
    os.makedirs(path_output)

# get image header information
phase_img = nb.load(phase_input)
phase_img.header["dim"][0] = 3
phase_img.header["dim"][4] = 1
header = phase_img.header
affine = phase_img.affine

phase_array = phase_img.get_fdata()

# normalize
phase_array = 2 * ( phase_array - np.min(phase_array) ) / ( np.max(phase_array) - np.min(phase_array) ) - 1
phase_array *= np.pi
phase_array = np.std(phase_array, axis=3)
phase_array[phase_array > phase_threshold] = 0
phase_array[phase_array != 0] = 1
phase_array -= 1
phase_array = np.abs(phase_array)

output = nb.Nifti1Image(phase_array, affine, header)
fileOUT = os.path.join(path_output,"vein.nii")
nb.save(output,fileOUT)

# write log
fileID = open(os.path.join(path_output,"vein_info.txt"),"a")
fileID.write("script executed: "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
fileID.write("phase_threshold: "+str(phase_threshold)+"\n")
fileID.close()
