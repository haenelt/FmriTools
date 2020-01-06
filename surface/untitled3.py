"""
Apply nonlinear deformation of surface to target space

The purpose of the following script is to deform a list or surface meshs to a target space using
a coordinate mapping.

Before running the script, login to queen via ssh and set the freesurfer environment by calling 
FREESURFER in the terminal.

created by Daniel Haenelt
Date created: 07-02-2019
Last modified: 09-12-2019
"""
from os.path import join, basename, splitext
from lib.io import get_filename
from lib.cmap import remove_edge_cmap
from lib.surface import deform_surface, remove_vertex_outliers
from lib.io.get_filename import get_filename
from lib.registration.get_scanner_transform import get_scanner_transform

# input files
input_surf = "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/lh.layer10"
input_orig = "/data/pt_01880/Experiment1_ODC/p3/anatomy/freesurfer/mri/orig.mgz"
input_ana = "/data/pt_01880/Experiment1_ODC/p3/anatomy/S23_MP2RAGE_0p7_UNI_Images_2.45.nii"
input_target = "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/diagnosis/mean_data.nii"
input_deform1 = "/data/pt_01880/odc_temp/deformation/header/T1_2_orig_scanner.nii" # ana2orig
input_deform2 = "/data/pt_01880/odc_temp/deformation/ge_epi2/epi2ana.nii.gz" # epi2ana
path_output = "/data/pt_01880/odc_temp/"

apply_header = True
remove_edge = False
remove_outliers = False

# to do
# probably trilinear interpolation only makes it nasty (do not do this)
# delete edges?
# get right output names
# deformation as iterations with all parameters as arrays
# write unsmoothed surface separately

""" do not edit below """

# get hemisphere from input surface file name
hemi = splitext(basename(input_surf))[0]

#%%

# deform surface
deform_surface(input_surf, 
               input_orig, 
               input_deform1, 
               input_ana,
               hemi, 
               path_output, 
               "trilinear",
               0,
               False,
               True)

#%%

deform_surface("/data/pt_01880/odc_temp/lh.layer10_def", 
               input_ana,
               input_deform2, 
               input_target,
               hemi,
               path_output, 
               "nearest",
               0,
               True,
               False)