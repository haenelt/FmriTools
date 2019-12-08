"""
Nonlinear deformation of surface to target space

The purpose of the following script is to deform a list or surface meshs to a target space using
a coordinate mapping.

Before running the script, login to queen via ssh and set the freesurfer environment by calling 
FREESURFER in the terminal.

created by Daniel Haenelt
Date created: 07-02-2019
Last modified: 08-12-2019
"""
import os
from lib.surface import deform_surface, remove_vertex_outliers

# input files
input_surf = ["/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/lh.layer0",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/lh.layer1",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/lh.layer2",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/lh.layer3",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/lh.layer4",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/lh.layer5",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/lh.layer6",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/lh.layer7",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/lh.layer8",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/lh.layer9",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/lh.layer10",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/rh.layer0",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/rh.layer1",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/rh.layer2",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/rh.layer3",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/rh.layer4",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/rh.layer5",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/rh.layer6",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/rh.layer7",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/rh.layer8",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/rh.layer9",
              "/data/pt_01880/Experiment1_ODC/p3/anatomy/layer/rh.layer10",
              ]

input_orig = "/data/pt_01880/Experiment1_ODC/p3/anatomy/freesurfer/mri/orig.mgz"
input_deform = "/data/pt_01880/Experiment1_ODC/p3/deformation/odc/GE_EPI2_rigid/epi2orig.nii.gz" # epi2orig
input_target = "/data/pt_01880/Experiment1_ODC/p3/odc/GE_EPI2/diagnosis/mean_data.nii"
path_output = "/data/pt_01880/Experiment1_ODC/p3/deformation/odc/GE_EPI2_rigid/epi_surf"

""" do not edit below """

for i in range(len(input_surf)):
    
    # get hemisphere from input surface file name
    hemi = os.path.splitext(os.path.basename(input_surf[i]))[0]
    
    # deform surface
    deform_surface(input_surf[i], 
                   input_orig, 
                   input_deform, 
                   input_target, 
                   hemi, 
                   path_output, 
                   False)

    # remove vertex outliers
    remove_vertex_outliers(os.path.join(path_output,os.path.basename(input_surf)+"_def"), 
                           os.path.join(path_output,os.path.basename(input_surf)+"_ind.txt"), 
                           5, 
                           True)