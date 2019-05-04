"""
FLASH -> orig registration

The purpose of the following script is to compute the deformation field for the registration 
between a GRE slab and the freesurfer orig. The script consists of the following
steps:
    1. set output folder structure
    2. scanner transform flash <-> mp2rage, mp2rage -> orig
    3. generate flash coordinate mapping
    4. apply scanner transform: mp2rage -> flash
    5. flirt registration: flash -> mp2rage
    6. apply flirt to flash coordinate mapping
    7. apply scanner transform to flash coordinate mapping
    8. apply final deformation to flash

Before running the script, login to queen via ssh and set the FSL environment by calling FSL in the 
terminal.

created by Daniel Haenelt
Date created: 18-04-2019
Last modified: 04-05-2019
"""
import os
import shutil as sh
from nipype.interfaces.fsl import FLIRT
from nipype.interfaces.fsl.preprocess import ApplyXFM
from nipype.interfaces.freesurfer.preprocess import MRIConvert
from nighres.registration import apply_coordinate_mappings
from lib.cmap import generate_coordinate_mapping
from lib.registration.get_scanner_transform import get_scanner_transform

# input data
file_flash = "/data/pt_01880/V2STRIPES/p6/anatomy/flash/S5_3D_GRE_3ech_iso0p5_slab_25.9.nii"
file_uni = "/data/pt_01880/V2STRIPES/p6/anatomy/UNI_0p7.nii"
file_orig = "/data/pt_01880/V2STRIPES/p6/anatomy/freesurfer/mri/orig.mgz"
path_output = "/data/pt_01880/V2STRIPES/p6/deformation/flash"
cleanup = False

""" do not edit below """

"""
set folder structure
"""
path_temp = os.path.join(path_output,"temp")

if not os.path.exists(path_output):
    os.makedirs(path_output)

if not os.path.exists(path_temp):
    os.makedirs(path_temp)

# copy input files
sh.copyfile(file_uni, os.path.join(path_temp,"uni.nii"))
sh.copyfile(file_flash, os.path.join(path_temp,"flash.nii"))

"""
convert orig to nifti
"""
mc = MRIConvert()
mc.inputs.in_file = file_orig
mc.inputs.out_file = os.path.join(path_temp,"orig.nii")
mc.inputs.in_type = "mgz"
mc.inputs.out_type = "nii"
mc.run()

"""
scanner transformation
"""
get_scanner_transform(os.path.join(path_temp,"uni.nii"), os.path.join(path_temp,"flash.nii"), path_temp)
get_scanner_transform(os.path.join(path_temp,"flash.nii"), os.path.join(path_temp,"uni.nii"), path_temp)
get_scanner_transform(os.path.join(path_temp,"uni.nii"), os.path.join(path_temp,"orig.nii"), path_temp)

"""
generate coordinate mapping
"""
generate_coordinate_mapping(os.path.join(path_temp,"flash.nii"), 0, path_temp, "flash", False)

"""
scanner transform uni to flash
"""
apply_coordinate_mappings(os.path.join(path_temp,"uni.nii"), # input 
                          os.path.join(path_temp,"uni_2_flash_scanner.nii"), # cmap
                          interpolation = "linear", # nearest or linear
                          padding = "zero", # closest, zero or max
                          save_data = True, # save output data to file (boolean)
                          overwrite = True, # overwrite existing results (boolean)
                          output_dir = path_temp, # output directory
                          file_name = "uni_apply_scanner" # base name with file extension for output
                          )

"""
flirt flash to uni
"""
os.chdir(path_temp)
flirt = FLIRT()
flirt.inputs.cost_func = "mutualinfo"
flirt.inputs.dof = 6
flirt.inputs.interp = "trilinear" # trilinear, nearestneighbour, sinc or spline
flirt.inputs.in_file = os.path.join(path_temp,"flash.nii")
flirt.inputs.reference = os.path.join(path_temp,"uni_apply_scanner_def-img.nii.gz")
flirt.inputs.output_type = "NIFTI_GZ"
flirt.inputs.out_file = os.path.join(path_temp,"flash_apply_flirt_def-img.nii.gz")
flirt.inputs.out_matrix_file = os.path.join(path_temp,"flirt_matrix.txt")
flirt.run()

"""
apply flirt to flash cmap
"""
applyxfm = ApplyXFM()
applyxfm.inputs.in_file = os.path.join(path_temp,"cmap_flash.nii")
applyxfm.inputs.reference = os.path.join(path_temp,"flash.nii")
applyxfm.inputs.in_matrix_file = os.path.join(path_temp,"flirt_matrix.txt")
applyxfm.inputs.interp = "trilinear"
applyxfm.inputs.padding_size = 0
applyxfm.inputs.output_type = "NIFTI_GZ"
applyxfm.inputs.out_file = os.path.join(path_temp,"cmap_apply_flirt_def-img.nii.gz")
applyxfm.inputs.apply_xfm = True
applyxfm.run() 

"""
apply scanner transform to flash cmap
"""
apply_coordinate_mappings(os.path.join(path_temp,"cmap_apply_flirt_def-img.nii.gz"),
                          os.path.join(path_temp,"flash_2_uni_scanner.nii"), # cmap 1
                          os.path.join(path_temp,"uni_2_orig_scanner.nii"), # cmap 2
                          interpolation = "linear", # nearest or linear
                          padding = "zero", # closest, zero or max
                          save_data = True, # save output data to file (boolean)
                          overwrite = True, # overwrite existing results (boolean)
                          output_dir = path_temp, # output directory
                          file_name = "cmap_apply_scanner" # base name with file extension for output
                          )

"""
apply deformation to source image
"""
apply_coordinate_mappings(os.path.join(path_temp,"flash.nii"), # input 
                          os.path.join(path_temp,"cmap_apply_scanner_def-img.nii.gz"),
                          interpolation = "linear", # nearest or linear
                          padding = "zero", # closest, zero or max
                          save_data = True, # save output data to file (boolean)
                          overwrite = True, # overwrite existing results (boolean)
                          output_dir = path_temp, # output directory
                          file_name = "flash_apply_deformation" # base name with file extension for output
                          )

# rename final deformation examples
os.rename(os.path.join(path_temp,"cmap_apply_scanner_def-img.nii.gz"),
          os.path.join(path_output,"flash2orig.nii.gz"))
os.rename(os.path.join(path_temp,"flash_apply_deformation_def-img.nii.gz"),
          os.path.join(path_output,"flash2orig_example.nii.gz"))

# clean intermediate files
if cleanup:
    sh.rmtree(path_temp, ignore_errors=True)
