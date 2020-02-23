"""
Vessel distance map

This scripts calculates a distance map to vessels found as ridge structures in a T2s weighted GRE.

created by Daniel Haenelt
Date created: 23-02-2020 
Last modified: 23-02-2020  
"""
import os
import numpy as np
import nibabel as nb
from nighres.filtering import filter_ridge_structures
from nighres.surface import probability_to_levelset

# list of GRE echoes
file_in = ["/home/daniel/Schreibtisch/temp/S13_3D_GRE_3ech_iso0p5_slab_8.42.nii",
           "/home/daniel/Schreibtisch/temp/S13_3D_GRE_3ech_iso0p5_slab_16.03.nii",
           "/home/daniel/Schreibtisch/temp/S13_3D_GRE_3ech_iso0p5_slab_25.nii",
           ]

# parameters
path_output = "/home/daniel/Schreibtisch"
img_res = 0.5
vessel_threshold = 0.3

""" do not edit below """

# initialize array
flash = nb.load(file_in[0])
flash_array = np.zeros_like(flash.get_fdata())

# get geometry average
for i in range(len(file_in)):
    flash_array += nb.load(file_in[i]).get_fdata() ** 2

flash_array = np.sqrt(flash_array)
flash = nb.Nifti1Image(flash_array, flash.affine, flash.header)

# filter ridge structures
ridge = filter_ridge_structures(flash,
                                structure_intensity='dark', 
                                output_type='probability', 
                                use_strict_min_max_filter=True, 
                                save_data=False, 
                                overwrite=False, 
                                output_dir=None, 
                                file_name=None)

# get binary vessel mask
ridge_array = ridge["result"].get_fdata()
ridge_array[ridge_array < vessel_threshold] = 0
ridge_array[ridge_array != 0] = 1
ridge = nb.Nifti1Image(ridge_array, ridge["result"].affine, ridge["result"].header)

# get from each pixel the distance to the closest vein pixel
vessel_distance = probability_to_levelset(ridge, 
                                          mask_image=None, 
                                          save_data=False,
                                          overwrite=False, 
                                          output_dir=None, 
                                          file_name=None)

# transform distance map to mm
vessel_distance_array = vessel_distance["result"].get_fdata()
vessel_distance_array *= 0.5 # vox2mm
vessel_distance_array[vessel_distance_array < 0] = 0 # remove negative distances
vessel_distance = nb.Nifti1Image(vessel_distance_array, 
                                 vessel_distance["result"].affine, 
                                 vessel_distance["result"].header)

# write output
nb.save(flash, os.path.join(path_output, "magn_geom_average.nii"))
nb.save(ridge, os.path.join(path_output,"vessel_filter.nii"))
nb.save(vessel_distance, os.path.join(path_output,"vessel_distance.nii"))