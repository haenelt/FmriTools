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

# input files
input_surf = ["/home/daniel/Schreibtisch/temo/lh.layer0",
              ]

input_orig = "/home/daniel/Schreibtisch/temo/orig.mgz"
input_deform = "/home/daniel/Schreibtisch/temo/epi2orig_edge.nii.gz" # epi2orig
input_target = "/home/daniel/Schreibtisch/temo/mean_data.nii"
path_output = "/home/daniel/Schreibtisch/output"

remove_edge = False

""" do not edit below """

for i in range(len(input_surf)):
    
    # get hemisphere from input surface file name
    hemi = splitext(basename(input_surf[i]))[0]
    
    # get path, basename and file extension from coordionate map
    path_deform, name_deform, ext_deform = get_filename(input_deform)
    
    # remove edges from cmap
    if remove_edge:
        remove_edge_cmap(input_deform, 5, 5)
        
        # rename coordinate map basename
        name_deform += "_edge"
    
    # deform surface
    deform_surface(input_surf[i], 
                   input_orig, 
                   join(path_deform, name_deform+ext_deform), 
                   input_target, 
                   hemi, 
                   path_output, 
                   True)
    
    # remove vertex outliers
    remove_vertex_outliers(join(path_output,basename(input_surf[i])+"_def"), 
                           join(path_output,basename(input_surf[i])+"_ind.txt"), 
                           5, 
                           True)
