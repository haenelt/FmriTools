# -*- coding: utf-8 -*-

# python standard library inputs
import os
import shutil as sh

# external inputs
import numpy as np
import nibabel as nb
from nighres.registration import apply_coordinate_mappings
from nighres.surface import probability_to_levelset
from nipype.interfaces import afni

# local inputs
from fmri_tools.processing import estimate_pv


"""
Partial volume estimation

This scripts calculates the partial volume contribution of wm, gm and csf based 
on computed wm (wm: 1, rest: 0) and csf (wm+gm: 1, rest: 0) masks. Optionally, 
the masks can be deformed to an target epi image. Levelset images are computed 
and upsampled. Partial voluming in the target epi space is computed by counting 
the tissue compartments in each epi voxel using a moving-average like algorithm.

The script needs an installation of afni.

created by Daniel Haenelt
Date created: 02-05-2019             
Last modified: 13-10-2020
"""

# input
file_wm = "/data/pt_01880/Experiment1_ODC/p4/anatomy/layer/left/wm_label.nii"
file_csf = "/data/pt_01880/Experiment1_ODC/p4/anatomy/layer/left/csf_label.nii"
file_target = "/data/pt_01880/Experiment1_ODC/p4/resting_state3/mean_udata.nii"
file_deformation = None
path_output = "/data/pt_01880/test_pve"
cleanup = False

# upsampled voxel size in mm
x_upsample = 0.2
y_upsample = 0.2
z_upsample = 0.2

# do not edit below

# make folder structure
path_levelset = os.path.join(path_output,"levelset")
path_border = os.path.join(path_output,"border")

if not os.path.exists(path_output):
    os.makedirs(path_output)

if not os.path.exists(path_levelset):
    os.makedirs(path_levelset)
    
if not os.path.exists(path_border):
    os.makedirs(path_border)

# deformation
if file_deformation is not None:
    wm_label = apply_coordinate_mappings(file_wm, # input 
                                         file_deformation, # cmap
                                         interpolation = "nearest", # nearest or linear
                                         padding = "zero", # closest, zero or max
                                         save_data = False, # save output data to file (boolean)
                                         overwrite = False, # overwrite existing results (boolean)
                                         output_dir = None, # output directory
                                         file_name = None # base name with file extension for output
                                         )
    
    csf_label = apply_coordinate_mappings(file_csf, # input 
                                          file_deformation, # cmap
                                          interpolation = "nearest", # nearest or linear
                                          padding = "zero", # closest, zero or max
                                          save_data = False, # save output data to file (boolean)
                                          overwrite = False, # overwrite existing results (boolean)
                                          output_dir = None, # output directory
                                          file_name = None # base name with file extension for output
                                          )
    
    wm_label = wm_label["result"]
    csf_label = csf_label["result"]
else:
    wm_label = nb.load(file_wm)
    csf_label = nb.load(file_csf)
    
# probability to levelset
probability_to_levelset(wm_label, 
                        save_data=True, 
                        overwrite=True, 
                        output_dir=path_levelset, 
                        file_name="wm_border")

probability_to_levelset(csf_label, 
                        save_data=True, 
                        overwrite=True, 
                        output_dir=path_levelset, 
                        file_name="csf_border")

# upsample levelsets
resample = afni.Resample()
resample.inputs.in_file = os.path.join(path_levelset,"wm_border_p2l-surf.nii.gz")
resample.inputs.outputtype = "NIFTI_GZ"
resample.inputs.out_file = os.path.join(path_levelset,"wm_border_p2l-surf_upsampled.nii.gz")
resample.inputs.resample_mode = "Li"
resample.inputs.voxel_size = (x_upsample,y_upsample,z_upsample)
resample.run()  

resample = afni.Resample()
resample.inputs.in_file = os.path.join(path_levelset,"csf_border_p2l-surf.nii.gz")
resample.inputs.outputtype = "NIFTI_GZ"
resample.inputs.out_file = os.path.join(path_levelset,"csf_border_p2l-surf_upsampled.nii.gz")
resample.inputs.resample_mode = "Li"
resample.inputs.voxel_size = (x_upsample,y_upsample,z_upsample)
resample.run()  

# binary border for gm, wm and csf
wm = nb.load(os.path.join(path_levelset,"wm_border_p2l-surf_upsampled.nii.gz"))
wm_array = wm.get_fdata()

wm_array[wm_array==0] = 1
wm_array[wm_array>0] = 0
wm_array[wm_array<0] = 1

csf = nb.load(os.path.join(path_levelset,"csf_border_p2l-surf_upsampled.nii.gz"))
csf_array = csf.get_fdata()

csf_array[csf_array>=0] = 1
csf_array[csf_array<0] = 0

gm_array = 1 - (wm_array + csf_array)

# After thresholding to binary images, CSF and WM have several common voxels. These are identified 
# in the GM image having voxel value -1. These voxels are changed to NaN in all binary images.
# Values GM = -1 are those voxels, which overlap in WM and CSF. These voxels are set to NaN.
wm_array[gm_array == -1] = np.nan
csf_array[gm_array == -1] = np.nan
gm_array[gm_array == -1] = np.nan

output = nb.Nifti1Image(wm_array,wm.affine,wm.header)
nb.save(output,os.path.join(path_border,"wm.nii.gz"))

output = nb.Nifti1Image(csf_array,csf.affine,csf.header)
nb.save(output,os.path.join(path_border,"csf.nii.gz"))

output = nb.Nifti1Image(gm_array,csf.affine,csf.header)
nb.save(output,os.path.join(path_border,"gm.nii.gz"))

# pv estimation
name_input = ["wm", "csf", "gm"]
for i in range(len(name_input)):
    estimate_pv(file_target, 
                os.path.join(path_border,name_input[i]+".nii.gz"), 
                path_output, 
                name_input[i])

# clean intermediate files
if cleanup:
    sh.rmtree(path_levelset, ignore_errors=True)
    sh.rmtree(path_border, ignore_errors=True)
