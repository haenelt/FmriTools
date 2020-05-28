"""
Partial volume estimation

This scripts calculates the partial volume contribution of wm, gm and csf based on the cortical 
ribbon mask from the freesurfer segmentation. Optionally, WM and CSF masks can be used as input
directly if no ribbon file is given. First, the ribbon mask is deformed to native epi space using 
nearest neighbor interpolation. Binary masks for the GM/WM and GM/CSF boundary are created 
(right WM: 41, right GM: 42, left WM: 2, right WM: 3) and levelset images are computed from those 
masks (positive outwards and negative inwards). The levelset images are upsampled and the partial
voluming in the original epi space is computed by counting the tissue compartments in each epi
voxel using a moving-average like algorithm.

Before running the script, login to queen via ssh and set the afni environment by calling AFNI in 
the terminal.

created by Daniel Haenelt
Date created: 02-05-2019             
Last modified: 28-05-2019  
"""
import os
import shutil as sh
import numpy as np
import nibabel as nb
from nighres.registration import apply_coordinate_mappings
from nighres.surface import probability_to_levelset
from nipype.interfaces import afni
from lib.processing import estimate_pv

# input
file_wm = None
file_csf = None
file_ribbon = "/data/pt_01880/V2STRIPES/p6/anatomy/freesurfer/mri/ribbon.mgz" # can be used instead of wm and csf
file_deformation = "/data/pt_01880/V2STRIPES/p6/deformation/resting_state/orig2epi.nii.gz"
path_output = "/data/pt_01880/V2STRIPES/p6/pve"
cleanup = False

# upsampled voxel size in mm
x_upsample = 0.2
y_upsample = 0.2
z_upsample = 0.2

""" do not edit below """

# make folder structure
path_ribbon = os.path.join(path_output,"ribbon")
path_binary = os.path.join(path_output,"binary")
path_levelset = os.path.join(path_output,"levelset")
path_border = os.path.join(path_output,"border")

if not os.path.exists(path_output):
    os.makedirs(path_output)

if not os.path.exists(path_ribbon):
    os.makedirs(path_ribbon)

if not os.path.exists(path_binary):
    os.makedirs(path_binary)
    
if not os.path.exists(path_levelset):
    os.makedirs(path_levelset)
    
if not os.path.exists(path_border):
    os.makedirs(path_border)

if file_ribbon is not None:
    
    # deform ribbon
    ribbon = apply_coordinate_mappings(file_ribbon, # input 
                                       file_deformation, # cmap
                                       interpolation = "nearest", # nearest or linear
                                       padding = "zero", # closest, zero or max
                                       save_data = False, # save output data to file (boolean)
                                       overwrite = False, # overwrite existing results (boolean)
                                       output_dir = None, # output directory
                                       file_name = None # base name with file extension for output
                                       )
    
    # create masks for binary white matter and csf borders
    data_array = ribbon["result"].get_fdata()
    
    wm_border = data_array.copy()
    wm_border[wm_border == 41] = 1
    wm_border[wm_border == 42] = 0
    wm_border[wm_border == 2] = 1
    wm_border[wm_border == 3] = 0
    
    output = nb.Nifti1Image(wm_border, ribbon["result"].affine, ribbon["result"].header)
    nb.save(output, os.path.join(path_binary,"wm_border.nii"))
    
    csf_border = data_array.copy()
    csf_border[csf_border == 41] = 1
    csf_border[csf_border == 42] = 1
    csf_border[csf_border == 2] = 1
    csf_border[csf_border == 3] = 1
    
    output = nb.Nifti1Image(csf_border, ribbon["result"].affine, ribbon["result"].header)
    nb.save(output, os.path.join(path_binary,"csf_border.nii"))
    
else:
    
    # deform wm
    wm_border = apply_coordinate_mappings(file_wm, # input 
                                          file_deformation, # cmap
                                          interpolation = "nearest", # nearest or linear
                                          padding = "zero", # closest, zero or max
                                          save_data = False, # save output data to file (boolean)
                                          overwrite = False, # overwrite existing results (boolean)
                                          output_dir = None, # output directory
                                          file_name = None # base name with file extension for output
                                          )
    
    nb.save(wm_border, os.path.join(path_binary,"wm_border.nii"))
    
    # deform csf
    csf_border = apply_coordinate_mappings(file_csf, # input 
                                           file_deformation, # cmap
                                           interpolation = "nearest", # nearest or linear
                                           padding = "zero", # closest, zero or max
                                           save_data = False, # save output data to file (boolean)
                                           overwrite = False, # overwrite existing results (boolean)
                                           output_dir = None, # output directory
                                           file_name = None # base name with file extension for output
                                           )
    
    nb.save(csf_border, os.path.join(path_binary,"csf_border.nii"))

# probability to levelset
probability_to_levelset(os.path.join(path_binary,"wm_border.nii"), 
                        save_data=True, 
                        overwrite=True, 
                        output_dir=path_levelset, 
                        file_name="wm_border")

probability_to_levelset(os.path.join(path_binary,"csf_border.nii"), 
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
    estimate_pv(os.path.join(path_ribbon,"ribbon_def-img.nii.gz"), 
                os.path.join(path_border,name_input[i]+".nii.gz"), 
                path_output, 
                name_input[i])

# clean intermediate files
if cleanup:
    sh.rmtree(path_ribbon, ignore_errors=True)
    sh.rmtree(path_binary, ignore_errors=True)
    sh.rmtree(path_levelset, ignore_errors=True)
    sh.rmtree(path_border, ignore_errors=True)