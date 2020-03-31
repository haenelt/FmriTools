"""
Map2ana

In the following script, source data is transformed using a deformation field and the transformed
data is mapped onto a surface mesh in the target space. The deformation of the input data is 
optional.

Before running the script, login to queen via ssh and set the freesurfer environments by calling 
FREESURFER in the terminal.

created by Daniel Haenelt
Date created: 06-03-2019            
Last modified: 31-01-2020  
"""
import os
from nighres.registration import apply_coordinate_mappings
from lib.io.get_filename import get_filename
from lib.mapping import map2surface

# input
file_in = ["/data/pt_01983/func/flicker/results/psc/native/psc_on_diff_GE_EPI2_diff.nii",
           "/data/pt_01983/func/flicker/results/psc/native/psc_on_GE_EPI2.nii",
           "/data/pt_01983/func/flicker/results/psc/native/psc_on_GE_EPI2_diff.nii",
           "/data/pt_01983/func/flicker/results/psc/native/psc_on_GE_EPI2_diff_calhoun.nii",
           "/data/pt_01983/func/flicker/results/psc/native/psc_on_GE_EPI2_tdm_early.nii",
           "/data/pt_01983/func/flicker/results/psc/native/psc_on_GE_EPI2_tdm_late.nii",
           "/data/pt_01983/func/flicker/results/psc/native/psc_on_off_diff_GE_EPI2_diff.nii",
           "/data/pt_01983/func/flicker/results/psc/native/psc_on_off_GE_EPI2.nii",
           "/data/pt_01983/func/flicker/results/psc/native/psc_on_off_GE_EPI2_diff.nii",
           "/data/pt_01983/func/flicker/results/psc/native/psc_on_off_GE_EPI2_diff_calhoun.nii",
           ]

surf_in = ["/data/pt_01983/func/anatomy/layer/lh.layer10",
           ]

# parameters
deformation_in1 = "/data/pt_01983/func/deformation/flicker/GE_EPI2/epi2orig.nii.gz"
deformation_in2 = None
path_output = "/data/pt_01983/func/flicker/results/psc"
interpolation = "linear" # can be linear or nearest

""" do not edit below """

# output folders
path_def = os.path.join(path_output,"def")
path_surf = os.path.join(path_output,"surf")

for i in range(len(file_in)):
    
    # apply deformation    
    if deformation_in1 is not None:
        apply_coordinate_mappings(file_in[i], 
                                  deformation_in1,
                                  deformation_in2, 
                                  interpolation=interpolation, 
                                  padding='closest', 
                                  save_data=True, 
                                  overwrite=True, 
                                  output_dir=path_def,
                                  file_name=None,
                                  )
    
    # get filename of deformed file
    if deformation_in1 is not None:
        _, name_file, _ = get_filename(file_in[i])
        filename_def = os.path.join(path_def, name_file+"_def-img.nii") 
    else:
        filename_def = file_in[i]
        
    # map to ana
    for j in range(len(surf_in)):
        
        # hemisphere
        hemi = os.path.splitext(os.path.basename(surf_in[j]))[0]
    
        # sample on surface
        map2surface(surf_in[j], 
                    filename_def, 
                    hemi, 
                    path_surf, 
                    input_white=None, 
                    input_ind=None, 
                    cleanup=True)