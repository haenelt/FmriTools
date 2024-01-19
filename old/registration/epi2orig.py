# -*- coding: utf-8 -*-
"""
EPI <-> ORIG registration

The purpose of the following script is to compute the deformation field for the
registration between antomy in conformed freesurfer space and native EPI space.
Inputs are a mean epi in native space, a T1 map of an mp2rage acquisition, a
freesurfer orig file generated from the mp2rage flat image and a skullstrip
mask. The mask does not have to be a binary image. A threshold value is defined
in the script which binarises the file. The mask is transformed to mp2rage space
if already in freesurfer space via scanner coordinates. The script consists of
the following steps:
    1. set output folder structure
    2. prepare mask
    3. convert orig.mgz to orig.nii
    4. scanner transform orig <-> t1 and t1 <-> epi
    5. n4 correction epi
    6. clean ana (remove ceiling and normalise)
    7. mask t1 and epi
    8. antsreg
    9. merge deformations
    10. clean deformations
    11. expand deformations
    12. apply deformations

The script needs an installation of freesurfer and ants.

"""

import os
import shutil as sh

import nibabel as nb
from nighres.registration import apply_coordinate_mappings, embedded_antsreg
from sh import gunzip

from ..io.vol import mri_convert
from ..registration.clean_ana import clean_ana
from ..registration.cmap import clean_coordinate_mapping, expand_coordinate_mapping
from ..registration.transform import apply_header, scanner_transform
from ..segmentation.mask import mask_ana, mask_epi
from ..utils.bias import remove_bias_ants

# input data
file_mean_epi = "/nobackup/actinium2/haenelt/ForOthers/RetinotopyFakhereh4/26958.af/retinotopy/diagnosis/mean_data.nii"
file_t1 = "/nobackup/actinium2/haenelt/ForOthers/RetinotopyFakhereh4/26958.af/anatomy/S7_MP2RAGE_0p7_T1_Images_2.45.nii"
file_mask = "/nobackup/actinium2/haenelt/ForOthers/RetinotopyFakhereh4/26958.af/anatomy/freesurfer/mri/brain.finalsurfs.mgz"  # brain.finalsurfs
file_orig = "/nobackup/actinium2/haenelt/ForOthers/RetinotopyFakhereh4/26958.af/anatomy/freesurfer/mri/orig.mgz"
path_output = "/nobackup/actinium2/haenelt/ForOthers/RetinotopyFakhereh4/26958.af/deformation/retinotopy"
clean_cmap = True
expand_cmap = True
cleanup = True

# parameters for mask preparation
mask_threshold = 1  # lower threshold for binarisation

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
cost_function = "CrossCorrelation"
interpolation = "Linear"

# do not edit below

# set folder structure
path_temp = os.path.join(path_output, "temp")
path_orig = os.path.join(path_temp, "orig")
path_scanner = os.path.join(path_temp, "scanner")
path_epi = os.path.join(path_temp, "epi")
path_t1 = os.path.join(path_temp, "t1")
path_syn = os.path.join(path_temp, "syn")

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
sh.copyfile(file_orig, os.path.join(path_orig, "orig.mgz"))
sh.copyfile(file_mean_epi, os.path.join(path_epi, "epi.nii"))
sh.copyfile(file_t1, os.path.join(path_t1, "T1.nii"))

# mask preparation

# convert to nifti
if os.path.splitext(file_mask)[1] == ".mgh":
    sh.copyfile(file_mask, os.path.join(path_t1, "mask.mgh"))
    mri_convert(os.path.join(path_t1, "mask.mgh"), os.path.join(path_t1, "mask.nii"))
elif os.path.splitext(file_mask)[1] == ".mgz":
    sh.copyfile(file_mask, os.path.join(path_t1, "mask.mgz"))
    mri_convert(os.path.join(path_t1, "mask.mgz"), os.path.join(path_t1, "mask.nii"))
elif os.path.splitext(file_mask)[1] == ".gz":
    sh.copyfile(file_mask, os.path.join(path_t1, "mask.nii.gz"))
    gunzip(os.path.join(path_t1, "mask.nii.gz"))
elif os.path.splitext(file_mask)[1] == ".nii":
    sh.copyfile(file_mask, os.path.join(path_t1, "mask.nii"))
else:
    raise ValueError("File extension of mask could not be identified!")

# binarise and overwrite mask
mask = nb.load(os.path.join(path_t1, "mask.nii"))
mask_array = mask.get_fdata()
mask_array[mask_array <= mask_threshold] = 0
mask_array[mask_array != 0] = 1
output = nb.Nifti1Image(mask_array, mask.affine, mask.header)
nb.save(output, os.path.join(path_t1, "mask.nii"))

# transform to mp2rage space via scanner coordinates
apply_header(
    os.path.join(path_t1, "mask.nii"),
    os.path.join(path_t1, "T1.nii"),
    os.path.join(path_t1, "mask.nii"),
    interp_method="nearest",
)

# convert orig to nifti
mri_convert(os.path.join(path_orig, "orig.mgz"), os.path.join(path_orig, "orig.nii"))

# scanner transformation
scanner_transform(
    os.path.join(path_orig, "orig.nii"),
    os.path.join(path_t1, "T1.nii"),
    path_scanner,
    False,
)
scanner_transform(
    os.path.join(path_t1, "T1.nii"),
    os.path.join(path_orig, "orig.nii"),
    path_scanner,
    False,
)

# bias field correction to epi
remove_bias_ants(os.path.join(path_epi, "epi.nii"), os.path.join(path_epi, "bepi.nii"))

# clean ana
clean_ana(os.path.join(path_t1, "T1.nii"), 1000.0, 4095.0, overwrite=True)

# mask t1 and epi
mask_ana(
    os.path.join(path_t1, "T1.nii"),
    os.path.join(path_t1, "mask.nii"),
    background_bright=False,
)
mask_epi(
    os.path.join(path_epi, "bepi.nii"),
    os.path.join(path_t1, "pT1.nii"),
    os.path.join(path_t1, "mask.nii"),
    niter_mask,
    sigma_mask,
)

# syn
embedded_antsreg(
    os.path.join(path_t1, "pT1.nii"),  # source image
    os.path.join(path_epi, "pbepi.nii"),  # target image
    run_rigid,  # whether or not to run a rigid registration first
    rigid_iterations,  # number of iterations in the rigid step
    run_affine,  # whether or not to run an affine registration first
    affine_iterations,  # number of iterations in the affine step
    run_syn,  # whether or not to run a SyN registration
    coarse_iterations,  # number of iterations at the coarse level
    medium_iterations,  # number of iterations at the medium level
    fine_iterations,  # number of iterations at the fine level
    cost_function,  # CrossCorrelation or MutualInformation
    interpolation,  # interpolation for registration result (NeareastNeighbor or Linear)
    convergence=1e-6,  # threshold for convergence (can make algorithm very slow)
    ignore_affine=False,  # ignore the affine matrix information extracted from the image header
    ignore_header=False,  # ignore the orientation information and affine matrix information extracted from the image header
    save_data=True,  # save output data to file
    overwrite=True,  # overwrite existing results
    output_dir=path_syn,  # output directory
    file_name="syn",  # output basename
)

# merge deformations
# orig -> epi
apply_coordinate_mappings(
    os.path.join(path_scanner, "orig_2_T1_scanner.nii"),  # input
    os.path.join(path_syn, "syn_ants-map.nii.gz"),  # cmap
    interpolation="linear",  # nearest or linear
    padding="zero",  # closest, zero or max
    save_data=True,  # save output data to file (boolean)
    overwrite=True,  # overwrite existing results (boolean)
    output_dir=path_output,  # output directory
    file_name="orig2epi",  # base name with file extension for output
)

# epi -> orig
apply_coordinate_mappings(
    os.path.join(path_syn, "syn_ants-invmap.nii.gz"),  # input
    os.path.join(path_scanner, "T1_2_orig_scanner.nii"),  # cmap
    interpolation="linear",  # nearest or linear
    padding="zero",  # closest, zero or max
    save_data=True,  # save output data to file (boolean)
    overwrite=True,  # overwrite existing results (boolean)
    output_dir=path_output,  # output directory
    file_name="epi2orig",  # base name with file extension for output
)

# rename final deformations
os.rename(
    os.path.join(path_output, "orig2epi_def-img.nii.gz"),
    os.path.join(path_output, "orig2epi.nii.gz"),
)
os.rename(
    os.path.join(path_output, "epi2orig_def-img.nii.gz"),
    os.path.join(path_output, "epi2orig.nii.gz"),
)

# clean deformation
if clean_cmap:
    epi2ana_cleaned = clean_coordinate_mapping(
        os.path.join(path_output, "orig2epi.nii.gz"),
        os.path.join(path_output, "epi2orig.nii.gz"),
        overwrite_file=True,
        save_mask=False,
    )

    # write mask
    nb.save(epi2ana_cleaned["mask"], os.path.join(path_output, "epi2orig_mask.nii.gz"))

# expand deformation
if expand_cmap:
    _ = expand_coordinate_mapping(
        os.path.join(path_output, "orig2epi.nii.gz"),
        path_output,
        name_output="orig2epi",
        write_output=True,
    )

    _ = expand_coordinate_mapping(
        os.path.join(path_output, "epi2orig.nii.gz"),
        path_output,
        name_output="epi2orig",
        write_output=True,
    )

# apply deformation
# orig -> epi
apply_coordinate_mappings(
    file_orig,  # input
    os.path.join(path_output, "orig2epi.nii.gz"),  # cmap
    interpolation="linear",  # nearest or linear
    padding="zero",  # closest, zero or max
    save_data=True,  # save output data to file (boolean)
    overwrite=True,  # overwrite existing results (boolean)
    output_dir=path_output,  # output directory
    file_name="orig2epi_example",  # base name with file extension for output
)

# epi -> orig
apply_coordinate_mappings(
    file_mean_epi,  # input
    os.path.join(path_output, "epi2orig.nii.gz"),  # cmap
    interpolation="linear",  # nearest or linear
    padding="zero",  # closest, zero or max
    save_data=True,  # save output data to file (boolean)
    overwrite=True,  # overwrite existing results (boolean)
    output_dir=path_output,  # output directory
    file_name="epi2orig_example",  # base name with file extension for output
)

# rename final deformation examples
os.rename(
    os.path.join(path_output, "orig2epi_example_def-img.nii.gz"),
    os.path.join(path_output, "orig2epi_example.nii.gz"),
)
os.rename(
    os.path.join(path_output, "epi2orig_example_def-img.nii.gz"),
    os.path.join(path_output, "epi2orig_example.nii.gz"),
)

# clean intermediate files
if cleanup:
    sh.rmtree(path_temp, ignore_errors=True)
