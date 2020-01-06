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

from lib.io.mgh2nii import mgh2nii
from lib.registration.get_scanner_transform import get_scanner_transform

# input files
input_surf = ["/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/lh.layer10",
              ]

input_orig = "/data/pt_01880/Experiment1_ODC/p3/anatomy/freesurfer/mri/orig.mgz"
input_ana = "/data/pt_01880/Experiment1_ODC/p3/anatomy/S23_MP2RAGE_0p7_UNI_Images_2.45.nii"
input_deform = "/data/pt_01880/odc_temp/deformation/ge_epi2/epi2ana.nii.gz" # epi2orig
input_target = "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/diagnosis/mean_data.nii"
path_output = "/data/pt_01880/odc_temp/"

apply_header = True
remove_edge = False
remove_outliers = False

""" do not edit below """

for i in range(len(input_surf)):
    
    # get hemisphere from input surface file name
    hemi = splitext(basename(input_surf[i]))[0]
    
    # get path, basename and file extension from coordionate map
    path_deform, name_deform, ext_deform = get_filename(input_deform)
    
    # apply header transformation to mesh
    if apply_header:
        mgh2nii(input_orig, path_output)
        get_scanner_transform(input_ana,join(path_output,"orig.nii"),path_output)
    
    # remove edges from cmap
    if remove_edge:
        remove_edge_cmap(input_deform, 5, 5)
        
        # rename coordinate map basename
        name_deform += "_edge"
    
    # deform surface
    deform_surface(input_surf[i], 
                   input_orig, 
                   "/data/pt_01880/odc_temp/S23_MP2RAGE_0p7_UNI_Images_2.45_2_orig_scanner.nii",#join(path_deform, name_deform+ext_deform), 
                   input_target, 
                   hemi, 
                   path_output, 
                   True)
    
    # remove vertex outliers
    if remove_outliers:
        remove_vertex_outliers(join(path_output,basename(input_surf[i])+"_def"), 
                               join(path_output,basename(input_surf[i])+"_ind.txt"), 
                               5, 
                               True)