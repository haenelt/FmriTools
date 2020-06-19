"""
MPM <-> MP2RAGE registration

The purpose of the following script is to compute the deformation field for the registration 
between MPM and MP2RAGE T1-maps. The script consists of the following steps:
    1. make folder structure
    2. scale PD-map (MPM)
    3. skullstrip MP2RAGE and MPM
    4. scale T1-map (MPM)
    5. scanner transformation (MPM <-> MP2RAGE)
    6. rigid registration (MPM <-> MP2RAGE)
    7. merge both transformations
    8. expand both transformations
    9. apply final transformation to MPM and MP2RAGE

created by Daniel Haenelt
Date created: 31-01-2019
Last modified: 19-06-2020
"""
import os
import shutil as sh
import numpy as np
import nibabel as nb
from nighres.registration import embedded_antsreg, apply_coordinate_mappings
from lib.cmap.expand_coordinate_mapping import expand_coordinate_mapping
from lib.registration.get_scanner_transform import get_scanner_transform
from lib.skullstrip.skullstrip_spm12 import skullstrip_spm12

# input files
file_mpm_r1 = "/home/daniel/Schreibtisch/mpm_reg/data/mpm/s122880a-145639-00001-00352-1_R1.nii"
file_mpm_pd = "/home/daniel/Schreibtisch/mpm_reg/data/mpm/s122880a-145639-00001-00352-1_PD.nii"
file_mp2rage_t1 = "/home/daniel/Schreibtisch/mpm_reg/data/mp2rage/S5_mp2rage_whole_brain_T1_Images_2.45.nii"
file_mp2rage_pd = "/home/daniel/Schreibtisch/mpm_reg/data/mp2rage/S4_mp2rage_whole_brain_INV2_2.45.nii"
expand_cmap = True
cleanup = False

# output path
path_output = "/home/daniel/Schreibtisch"

# set environments
pathSPM12 = "/home/daniel/source/spm12"

# parameters for syn 
run_rigid = True
rigid_iterations = 1000 
run_affine = False 
affine_iterations = 1000 
run_syn = False
coarse_iterations = 50 
medium_iterations = 150 
fine_iterations = 100 
cost_function = 'CrossCorrelation' 
interpolation = 'Linear'

""" do not edit below """

# make folder structure
print("make folder structure")

path_res = os.path.join(path_output,"registration")
path_mpm = os.path.join(path_res,"mpm")
path_mp2rage = os.path.join(path_res,"mp2rage")
path_deformation = os.path.join(path_res,"deformation")

if not os.path.exists(path_res):
    os.makedirs(path_res)

if not os.path.exists(path_mpm):
    os.mkdir(path_mpm)

if not os.path.exists(path_mp2rage):
    os.mkdir(path_mp2rage)

if not os.path.exists(path_deformation):
    os.mkdir(path_deformation)

# prepare PD-map (MPM)
print("prepare mpm (pd)")

mpm_img = nb.load(file_mpm_pd)
mpm_array = mpm_img.get_fdata()
mpm_array[mpm_array <= 0] = 0
mpm_array = mpm_array.astype('uint16')

mpm_img.set_data_dtype('uint16')
output = nb.Nifti1Image(mpm_array, mpm_img.affine, mpm_img.header)
nb.save(output,os.path.join(path_mpm,"mpm_pd.nii"))

# skull stripping
print("skullstrip mp2rage")
skullstrip_spm12(file_mp2rage_pd, pathSPM12, path_mp2rage)

print("skullstrip mpm")
skullstrip_spm12(os.path.join(path_mpm,"mpm_pd.nii"), pathSPM12, path_mpm)

# prepare T1-map (MPM)
print("prepare mpm (t1)")

mask_img = nb.load(os.path.join(path_mpm,"skull","skullstrip_mask.nii"))
mask_array = mask_img.get_fdata()

mpm_img = nb.load(file_mpm_r1)
mpm_array = mpm_img.get_fdata()

mpm_array[mpm_array <= 0] = np.nan
mpm_array = 1000 / mpm_array
mpm_array[np.isnan(mpm_array)] = 0
mpm_array[mpm_array > 3999] = 3999
mpm_array = mpm_array * mask_array
mpm_array = mpm_array.astype('uint16')

mpm_img.set_data_dtype('uint16')
output = nb.Nifti1Image(mpm_array, mpm_img.affine, mpm_img.header)
nb.save(output,os.path.join(path_mpm,"mpm_t1.nii"))

# prepare T1-map (MP2RAGE)
print("prepare mp2rage (t1)")

mask_img = nb.load(os.path.join(path_mp2rage,"skull","skullstrip_mask.nii"))
mask_array = mask_img.get_fdata()

mp2rage_img = nb.load(file_mp2rage_t1)
mp2rage_array = mp2rage_img.get_fdata()

mp2rage_array = mp2rage_array * mask_array
mp2rage_array = mp2rage_array.astype('uint16')

mp2rage_img.set_data_dtype('uint16')
output = nb.Nifti1Image(mp2rage_array, mp2rage_img.affine, mp2rage_img.header)
nb.save(output,os.path.join(path_mp2rage,"mp2rage_t1.nii"))

# scanner transformation (MPM -> MP2RAGE)
print("get scanner transformation")
get_scanner_transform(os.path.join(path_mpm,"mpm_t1.nii"),
                      os.path.join(path_mp2rage,"mp2rage_t1.nii"), 
                      path_deformation, 
                      False)

get_scanner_transform(os.path.join(path_mp2rage,"mp2rage_t1.nii"), 
                      os.path.join(path_mpm,"mpm_t1.nii"),
                      path_deformation, 
                      False)

# apply scanner transformation to MPM
apply_coordinate_mappings(os.path.join(path_mpm,"mpm_t1.nii"), # input
                          os.path.join(path_deformation,"mpm_t1_2_mp2rage_t1_scanner.nii"), # first cmap
                          mapping2 = None, # second cmap
                          mapping3 = None, # third cmap
                          mapping4 = None, # fourth cmap
                          interpolation = "linear", # nearest or linear
                          padding = "zero", # closest, zero or max
                          save_data = True, # save output data to file (boolean)
                          overwrite = True, # overwrite existing results (boolean)
                          output_dir = path_mpm, # output directory
                          file_name = "scanner" # base name with file extension for output
                          )

# rigid registration (SYN)
embedded_antsreg(os.path.join(path_mpm,"scanner_def-img.nii.gz"), # source image
                 os.path.join(path_mp2rage,"mp2rage_t1.nii"), # target image 
                 run_rigid, # whether or not to run a rigid registration first 
                 rigid_iterations, # number of iterations in the rigid step
                 run_affine, # whether or not to run an affine registration first
                 affine_iterations, # number of iterations in the affine step
                 run_syn, # whether or not to run a SyN registration
                 coarse_iterations, # number of iterations at the coarse level
                 medium_iterations, # number of iterations at the medium level
                 fine_iterations, # number of iterations at the fine level
                 cost_function, # CrossCorrelation or MutualInformation
                 interpolation, # interpolation for registration result (NeareastNeighbor or Linear)
                 convergence = 1e-6, # threshold for convergence (can make algorithm very slow)
                 ignore_affine = False, # ignore the affine matrix information extracted from the image header 
                 ignore_header = False, # ignore the orientation information and affine matrix information extracted from the image header
                 save_data = True, # save output data to file
                 overwrite = True, # overwrite existing results 
                 output_dir = path_deformation, # output directory
                 file_name = "rigid", # output basename
                 )

# merge coordinate mappings
apply_coordinate_mappings(os.path.join(path_deformation,"mpm_t1_2_mp2rage_t1_scanner.nii"), # input 
                          os.path.join(path_deformation,"rigid_ants-map.nii.gz"), # first cmap
                          mapping2 = None, # second cmap
                          mapping3 = None, # third cmap
                          mapping4 = None, # fourth cmap
                          interpolation = "linear", # nearest or linear
                          padding = "zero", # closest, zero or max
                          save_data = True, # save output data to file (boolean)
                          overwrite = True, # overwrite existing results (boolean)
                          output_dir = path_output, # output directory
                          file_name = "mpm_2_mp2rage" # base name with file extension for output
                          )

apply_coordinate_mappings(os.path.join(path_deformation,"rigid_ants-invmap.nii.gz"), # input
                          os.path.join(path_deformation,"mp2rage_t1_2_mpm_t1_scanner.nii"), # first cmap
                          mapping2 = None, # second cmap
                          mapping3 = None, # third cmap
                          mapping4 = None, # fourth cmap
                          interpolation = "linear", # nearest or linear
                          padding = "zero", # closest, zero or max
                          save_data = True, # save output data to file (boolean)
                          overwrite = True, # overwrite existing results (boolean)
                          output_dir = path_output, # output directory
                          file_name = "mp2rage_2_mpm" # base name with file extension for output
                          )

# rename final deformation examples
os.rename(os.path.join(path_output,"mpm_2_mp2rage_def-img.nii.gz"),
          os.path.join(path_output,"mpm_2_mp2rage.nii.gz"))
os.rename(os.path.join(path_output,"mp2rage_2_mpm_def-img.nii.gz"),
          os.path.join(path_output,"mp2rage_2_mpm.nii.gz"))

"""
expand deformation
"""
if expand_cmap:
    _ = expand_coordinate_mapping(os.path.join(path_output, "mpm_2_mp2rage.nii.gz"),
                                  path_output, 
                                  name_output="mpm_2_mp2rage", 
                                  write_output=True)
    
    _ = expand_coordinate_mapping(os.path.join(path_output, "mp2rage_2_mpm.nii.gz"),
                                  path_output, 
                                  name_output="mp2rage_2_mpm", 
                                  write_output=True)

# apply final deformation to MPM
apply_coordinate_mappings(file_mpm_r1, # input 
                          os.path.join(path_output,"mpm_2_mp2rage.nii.gz"), # first cmap
                          mapping2 = None, # second cmap
                          mapping3 = None, # third cmap
                          mapping4 = None, # fourth cmap
                          interpolation = "linear", # nearest or linear
                          padding = "zero", # closest, zero or max
                          save_data = True, # save output data to file (boolean)
                          overwrite = True, # overwrite existing results (boolean)
                          output_dir = path_output, # output directory
                          file_name = "mpm_2_mp2rage_example" # base name with file extension for output
                          )

apply_coordinate_mappings(file_mp2rage_t1, # input 
                          os.path.join(path_output,"mp2rage_2_mpm.nii.gz"), # first cmap
                          mapping2 = None, # second cmap
                          mapping3 = None, # third cmap
                          mapping4 = None, # fourth cmap
                          interpolation = "linear", # nearest or linear
                          padding = "zero", # closest, zero or max
                          save_data = True, # save output data to file (boolean)
                          overwrite = True, # overwrite existing results (boolean)
                          output_dir = path_output, # output directory
                          file_name = "mp2rage_2_mpm_example" # base name with file extension for output
                          )

# rename final deformation examples
os.rename(os.path.join(path_output,"mpm_2_mp2rage_example_def-img.nii.gz"),
          os.path.join(path_output,"mpm_2_mp2rage_example.nii.gz"))
os.rename(os.path.join(path_output,"mp2rage_2_mpm_example_def-img.nii.gz"),
          os.path.join(path_output,"mp2rage_2_mpm_example.nii.gz"))

# clean intermediate files
if cleanup:
    sh.rmtree(path_res, ignore_errors=True)