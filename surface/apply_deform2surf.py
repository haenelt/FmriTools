"""
Apply nonlinear deformation to surface mesh

The purpose of the following script is to apply two deformations in succession (orig -> ana -> epi)
to a surface mesh using generated coordinate mappings.

Before running the script, login to queen via ssh and set the freesurfer environment by calling 
FREESURFER in the terminal.

created by Daniel Haenelt
Date created: 07-02-2019
Last modified: 07-01-2020
"""
import os 
from os.path import join, basename, splitext
from lib.surface import deform_surface

# input files
input_surf = ["/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/lh.layer10"]
input_orig = "/data/pt_01880/Experiment1_ODC/p3/anatomy/freesurfer/mri/orig.mgz"
input_ana = "/data/pt_01880/Experiment1_ODC/p3/anatomy/S23_MP2RAGE_0p7_UNI_Images_2.45.nii"
input_epi = "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/diagnosis/mean_data.nii"
input_deform1 = "/data/pt_01880/odc_temp/deformation/header/T1_2_orig_scanner.nii" # ana2orig
input_deform2 = "/data/pt_01880/odc_temp/deformation/ge_epi2/epi2ana.nii.gz" # epi2ana
path_output = "/data/pt_01880/odc_temp/"

""" do not edit below """

for i in range(len(input_surf)):    
    
    # orig -> ana
    deform_surface(input_surf[i], 
                   input_orig, 
                   input_deform1, 
                   input_ana, 
                   splitext(basename(input_surf[i]))[0], 
                   path_output, 
                   smooth_iter=10, 
                   sort_faces=False,
                   cleanup=True)

    # rename output
    os.rename(join(path_output, basename(input_surf[i])+"_def"),
              join(path_output, basename(input_surf[i])+"_def1"))
    os.rename(join(path_output, basename(input_surf[i])+"_def_smooth"),
              join(path_output, basename(input_surf[i])+"_def1_smooth")) 

    # ana -> epi
    deform_surface(join(path_output, basename(input_surf[i])+"_def1_smooth"),
                   input_ana,
                   input_deform2, 
                   input_epi,
                   splitext(basename(input_surf[i]))[0], 
                   path_output, 
                   smooth_iter=10, 
                   sort_faces=True,
                   cleanup=False)

    # rename output
    os.rename(join(path_output, basename(input_surf[i])+"_def1_smooth_def"),
              join(path_output, basename(input_surf[i])+"_def2"))
    os.rename(join(path_output, basename(input_surf[i])+"_def1_smooth_def_smooth"),
              join(path_output, basename(input_surf[i])+"_def2_smooth"))
    os.rename(join(path_output, basename(input_surf[i])+"_def1_smooth_ind.txt"),
              join(path_output, basename(input_surf[i])+"_def2_ind"))