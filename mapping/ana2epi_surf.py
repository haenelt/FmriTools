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
input_surf = ["/home/daniel/Schreibtisch/temo/lh.layer0",
              ]

input_orig = "/home/daniel/Schreibtisch/temo/orig.mgz"
input_deform = "/home/daniel/Schreibtisch/temo/epi2orig.nii.gz" # epi2orig
input_target = "/home/daniel/Schreibtisch/temo/mean_data.nii"
path_output = "/home/daniel/Schreibtisch/output"

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