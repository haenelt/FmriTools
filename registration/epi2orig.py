"""
EPI <-> ORIG registration

The purpose of the following script is to compute the deformation field for the registration 
between antomy in conformed freesurfer space and native EPI space. Inputs are a mean epi in native 
space, a T1 map of an mp2rage acquisition, a freesurfer orig file generated from the mp2rage flat
image and a skullstrip mask. The mask does not have to be a binary image. A threshold value is 
defined in the script which binarises the file. The mask is transformed to mp2rage space if already 
in freesurfer space via scanner coordinates. The script consists of the following steps:
    1. set output folder structure
    2. prepare mask
    3. convert orig.mgz to orig.nii
    4. scanner transform orig <-> t1 and t1 <-> epi
    5. n4 correction epi
    6. mask t1 and epi
    7. antsreg
    8. merge deformations
    9. apply deformations

Before running the script, login to queen via ssh and set the freesurfer and ANTS environments by 
calling FREESURFER and ANTSENV in the terminal.

created by Daniel Haenelt
Date created: 10-01-2019
Last modified: 15-07-2019
"""
import os
import shutil as sh
import nibabel as nb
from sh import gunzip
from nipype.interfaces.freesurfer import ApplyVolTransform
from nipype.interfaces.freesurfer.preprocess import MRIConvert
from nipype.interfaces.ants import N4BiasFieldCorrection
from nighres.registration import embedded_antsreg, apply_coordinate_mappings
from lib.registration.get_scanner_transform import get_scanner_transform
from lib.registration.mask_ana import mask_ana
from lib.registration.mask_epi import mask_epi

# input data
file_mean_epi = "/data/pt_01880/Experiment2_Rivalry/p3/retinotopy/diagnosis/mean_data.nii"
file_t1 = "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/S7_MP2RAGE_0p7_T1_Images_2.45.nii"
file_mask = "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/freesurfer/mri/brain.finalsurfs.manedit.mgz"
file_orig = "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/freesurfer/mri/orig.mgz"
path_output = "/data/pt_01880/Experiment2_Rivalry/p3/deformation/retinotopy"
cleanup = True

# parameters for mask preparation
mask_threshold = 1 # lower threshold for binarisation

# parameters for epi skullstrip
niter_mask = 3
sigma_mask = 3

# parameters for syn 
run_rigid = True
rigid_iterations = 1000 
run_affine = False 
affine_iterations = 1000 
run_syn = True 
coarse_iterations = 50 
medium_iterations = 150 
fine_iterations = 100 
cost_function = 'CrossCorrelation' 
interpolation = 'Linear' 

""" do not edit below """

"""
set folder structure
"""
path_temp = os.path.join(path_output,"temp")
path_orig = os.path.join(path_temp,"orig")
path_scanner = os.path.join(path_temp,"scanner")
path_epi = os.path.join(path_temp,"epi")
path_t1 = os.path.join(path_temp,"t1")
path_syn = os.path.join(path_temp,"syn")

if not os.path.exists(path_output):
    os.makedirs(path_output)

if not os.path.exists(path_temp):
    os.makedirs(path_temp)

if not os.path.exists(path_orig):
    os.makedirs(path_orig)

if not os.path.exists(path_scanner):
    os.makedirs(path_scanner)

if not os.path.exists(path_epi):
    os.makedirs(path_epi)

if not os.path.exists(path_t1):
    os.makedirs(path_t1)

if not os.path.exists(path_syn):
    os.makedirs(path_syn)

# copy input files
sh.copyfile(file_orig, os.path.join(path_orig,"orig.mgz"))
sh.copyfile(file_mean_epi, os.path.join(path_epi,"epi.nii"))
sh.copyfile(file_t1, os.path.join(path_t1,"T1.nii"))

"""
mask preparation
"""

# convert to nifti
if os.path.splitext(file_mask)[1] == ".mgh":
    sh.copyfile(file_mask, os.path.join(path_t1,"mask.mgh"))
    mc = MRIConvert()
    mc.inputs.in_file = os.path.join(path_t1,"mask.mgh")
    mc.inputs.out_file = os.path.join(path_t1,"mask.nii")
    mc.inputs.in_type = "mgh"
    mc.inputs.out_type = "nii"
    mc.run()    
elif os.path.splitext(file_mask)[1] == ".mgz":
    sh.copyfile(file_mask, os.path.join(path_t1,"mask.mgz"))
    mc = MRIConvert()
    mc.inputs.in_file = os.path.join(path_t1,"mask.mgz")
    mc.inputs.out_file = os.path.join(path_t1,"mask.nii")
    mc.inputs.in_type = "mgz"
    mc.inputs.out_type = "nii"
    mc.run()  
elif os.path.splitext(file_mask)[1] == ".gz":
    sh.copyfile(file_mask, os.path.join(path_t1,"mask.nii.gz"))
    gunzip(os.path.join(path_t1,"mask.nii.gz"))
elif os.path.splitext(file_mask)[1] == ".nii":
    sh.copyfile(file_mask, os.path.join(path_t1,"mask.nii"))
else:
    print("File extension of mask could not be identified!")

# binarise and overwrite mask
mask = nb.load(os.path.join(path_t1,"mask.nii"))
mask_array = mask.get_fdata()
mask_array[mask_array <= mask_threshold] = 0
mask_array[mask_array != 0] = 1
output = nb.Nifti1Image(mask_array, mask.affine, mask.header)
nb.save(output,os.path.join(path_t1,"mask.nii"))

# transform to mp2rage space via scanner coordinates
applyreg = ApplyVolTransform()
applyreg.inputs.source_file = os.path.join(path_t1,"mask.nii")
applyreg.inputs.target_file = os.path.join(path_t1,"T1.nii")
applyreg.inputs.transformed_file = os.path.join(path_t1,"mask.nii")
applyreg.inputs.reg_header = True
applyreg.inputs.interp = "nearest"
applyreg.run()

"""
convert orig to nifti
"""
mc = MRIConvert()
mc.inputs.in_file = os.path.join(path_orig,"orig.mgz")
mc.inputs.out_file = os.path.join(path_orig,"orig.nii")
mc.inputs.in_type = "mgz"
mc.inputs.out_type = "nii"
mc.run()

"""
scanner transformation
"""
get_scanner_transform(os.path.join(path_orig,"orig.nii"),os.path.join(path_t1,"T1.nii"),path_scanner)
get_scanner_transform(os.path.join(path_t1,"T1.nii"),os.path.join(path_orig,"orig.nii"),path_scanner)

"""
bias field correction to epi
"""
n4 = N4BiasFieldCorrection()
n4.inputs.dimension = 3
n4.inputs.input_image = os.path.join(path_epi,"epi.nii")
n4.inputs.bias_image = os.path.join(path_epi,'n4bias.nii')
n4.inputs.output_image = os.path.join(path_epi,"bepi.nii")
n4.run()

"""
mask t1 and epi
"""
mask_ana(os.path.join(path_t1,"T1.nii"),os.path.join(path_t1,"mask.nii"))
mask_epi(os.path.join(path_epi,"bepi.nii"), 
         os.path.join(path_t1,"pT1.nii"), 
         os.path.join(path_t1,"mask.nii"), 
         niter_mask, sigma_mask)

"""
syn
"""
embedded_antsreg(os.path.join(path_t1,"pT1.nii"), # source image
                 os.path.join(path_epi,"pbepi.nii"), # target image 
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
                 output_dir = path_syn, # output directory
                 file_name = "syn", # output basename
                 )

"""
merge deformations
"""
# orig -> epi
apply_coordinate_mappings(os.path.join(path_scanner,"orig_2_T1_scanner.nii"), # input 
                          os.path.join(path_syn,"syn_ants-map.nii.gz"), # cmap
                          interpolation = "linear", # nearest or linear
                          padding = "zero", # closest, zero or max
                          save_data = True, # save output data to file (boolean)
                          overwrite = True, # overwrite existing results (boolean)
                          output_dir = path_output, # output directory
                          file_name = "orig2epi" # base name with file extension for output
                          )

# epi -> orig
apply_coordinate_mappings(os.path.join(path_syn,"syn_ants-invmap.nii.gz"), # input
                          os.path.join(path_scanner,"T1_2_orig_scanner.nii"), # cmap
                          interpolation = "linear", # nearest or linear
                          padding = "zero", # closest, zero or max
                          save_data = True, # save output data to file (boolean)
                          overwrite = True, # overwrite existing results (boolean)
                          output_dir = path_output, # output directory
                          file_name = "epi2orig" # base name with file extension for output
                          )

# rename final deformations
os.rename(os.path.join(path_output,"orig2epi_def-img.nii.gz"),
          os.path.join(path_output,"orig2epi.nii.gz"))
os.rename(os.path.join(path_output,"epi2orig_def-img.nii.gz"),
          os.path.join(path_output,"epi2orig.nii.gz"))

"""
apply deformation
"""
# orig -> epi
apply_coordinate_mappings(file_orig, # input 
                          os.path.join(path_output,"orig2epi.nii.gz"), # cmap
                          interpolation = "linear", # nearest or linear
                          padding = "zero", # closest, zero or max
                          save_data = True, # save output data to file (boolean)
                          overwrite = True, # overwrite existing results (boolean)
                          output_dir = path_output, # output directory
                          file_name = "orig2epi_example" # base name with file extension for output
                          )

# epi -> orig
apply_coordinate_mappings(file_mean_epi, # input 
                          os.path.join(path_output,"epi2orig.nii.gz"), # cmap
                          interpolation = "linear", # nearest or linear
                          padding = "zero", # closest, zero or max
                          save_data = True, # save output data to file (boolean)
                          overwrite = True, # overwrite existing results (boolean)
                          output_dir = path_output, # output directory
                          file_name = "epi2orig_example" # base name with file extension for output
                          )

# rename final deformation examples
os.rename(os.path.join(path_output,"orig2epi_example_def-img.nii.gz"),
          os.path.join(path_output,"orig2epi_example.nii.gz"))
os.rename(os.path.join(path_output,"epi2orig_example_def-img.nii.gz"),
          os.path.join(path_output,"epi2orig_example.nii.gz"))

# clean intermediate files
if cleanup:
    sh.rmtree(path_temp, ignore_errors=True)
