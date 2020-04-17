"""
Apply deformation to surface mesh

The purpose of the following script is to apply two deformations in succession (orig -> ana -> mpm)
to a surface mesh using generated coordinate mappings.

Before running the script, login to queen via ssh and set the freesurfer environment by calling 
FREESURFER in the terminal.

created by Daniel Haenelt
Date created: 30-01-2020
Last modified: 10-03-2020
"""
import os 
from os.path import join, basename, splitext
from lib.surface import deform_surface

# input files
input_surf = ["/data/pt_01983/func/anatomy/layer/lh.layer0",
              "/data/pt_01983/func/anatomy/layer/lh.layer1",
              "/data/pt_01983/func/anatomy/layer/lh.layer2",
              "/data/pt_01983/func/anatomy/layer/lh.layer3",
              "/data/pt_01983/func/anatomy/layer/lh.layer4",
              "/data/pt_01983/func/anatomy/layer/lh.layer5",
              "/data/pt_01983/func/anatomy/layer/lh.layer6",
              "/data/pt_01983/func/anatomy/layer/lh.layer7",
              "/data/pt_01983/func/anatomy/layer/lh.layer8",
              "/data/pt_01983/func/anatomy/layer/lh.layer9",
              "/data/pt_01983/func/anatomy/layer/lh.layer10",
              "/data/pt_01983/func/anatomy/layer/lh.layer11",
              "/data/pt_01983/func/anatomy/layer/lh.layer12",
              "/data/pt_01983/func/anatomy/layer/lh.layer13",
              "/data/pt_01983/func/anatomy/layer/lh.layer14",
              "/data/pt_01983/func/anatomy/layer/lh.layer15",
              "/data/pt_01983/func/anatomy/layer/lh.layer16",
              "/data/pt_01983/func/anatomy/layer/lh.layer17",
              "/data/pt_01983/func/anatomy/layer/lh.layer18",
              "/data/pt_01983/func/anatomy/layer/lh.layer19",
              "/data/pt_01983/func/anatomy/layer/lh.layer20",
              "/data/pt_01983/func/anatomy/layer/rh.layer0",
              "/data/pt_01983/func/anatomy/layer/rh.layer1",
              "/data/pt_01983/func/anatomy/layer/rh.layer2",
              "/data/pt_01983/func/anatomy/layer/rh.layer3",
              "/data/pt_01983/func/anatomy/layer/rh.layer4",
              "/data/pt_01983/func/anatomy/layer/rh.layer5",
              "/data/pt_01983/func/anatomy/layer/rh.layer6",
              "/data/pt_01983/func/anatomy/layer/rh.layer7",
              "/data/pt_01983/func/anatomy/layer/rh.layer8",
              "/data/pt_01983/func/anatomy/layer/rh.layer9",
              "/data/pt_01983/func/anatomy/layer/rh.layer10",
              "/data/pt_01983/func/anatomy/layer/rh.layer11",
              "/data/pt_01983/func/anatomy/layer/rh.layer12",
              "/data/pt_01983/func/anatomy/layer/rh.layer13",
              "/data/pt_01983/func/anatomy/layer/rh.layer14",
              "/data/pt_01983/func/anatomy/layer/rh.layer15",
              "/data/pt_01983/func/anatomy/layer/rh.layer16",
              "/data/pt_01983/func/anatomy/layer/rh.layer17",
              "/data/pt_01983/func/anatomy/layer/rh.layer18",
              "/data/pt_01983/func/anatomy/layer/rh.layer19",
              "/data/pt_01983/func/anatomy/layer/rh.layer20",
              ]
input_orig = "/data/pt_01983/func/anatomy/freesurfer/mri/orig.mgz"
input_ana = "/data/pt_01983/func/anatomy/S7_MP2RAGE_0p7_T1_Images_2.45.nii"
input_mpm = "/data/pt_01983/func/mpm/Session1/seste/s2019-10-17_13-17-140642-00001-00352-1__dis3d_R1_scaled.nii"
input_deform1 = "/data/pt_01983/func/deformation/header/T12orig.nii.gz" # ana2orig
input_deform2 = "/data/pt_01983/func/deformation/MPM/sess1/mpm_2_mp2rage.nii.gz" # mpm2ana
path_output = "/data/pt_01983/func/anatomy/layer_mpm/Session1"

""" do not edit below """

for i in range(len(input_surf)):    
    
    # orig -> ana
    deform_surface(input_surf[i], 
                   input_orig, 
                   input_deform1, 
                   input_ana, 
                   splitext(basename(input_surf[i]))[0], 
                   path_output, 
                   interp_method="nearest",
                   smooth_iter=10, 
                   sort_faces=False,
                   flip_faces=False,
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
                   input_mpm,
                   splitext(basename(input_surf[i]))[0], 
                   path_output, 
                   interp_method="nearest",
                   smooth_iter=10, 
                   sort_faces=False,
                   flip_faces=False,
                   cleanup=True)

    # rename output
    os.rename(join(path_output, basename(input_surf[i])+"_def1_smooth_def"),
              join(path_output, basename(input_surf[i])+"_def2"))
    os.rename(join(path_output, basename(input_surf[i])+"_def1_smooth_def_smooth"),
              join(path_output, basename(input_surf[i])+"_def2_smooth"))
    
    # remove all intermediate steps
    os.remove(join(path_output, basename(input_surf[i])+"_def1"))
    os.remove(join(path_output, basename(input_surf[i])+"_def1_smooth"))
    os.remove(join(path_output, basename(input_surf[i])+"_def2"))

    # rename final surface
    os.rename(join(path_output, basename(input_surf[i])+"_def2_smooth"),
              join(path_output, basename(input_surf[i])+"_mpm"))
