"""
Map2ana

In the following script, source data is transformed using a deformation field and the transformed
data is mapped onto a surface mesh in the target space.

Before running the script, login to queen via ssh and set the freesurfer environments by calling 
FREESURFER in the terminal.

created by Daniel Haenelt
Date created: 06-03-2019            
Last modified: 02-09-2019  
"""
import os
from nighres.registration import apply_coordinate_mappings
from lib.mapping import map2surface

# input
input_file = [        
        "/data/pt_01880/Experiment2_Rivalry/p3/resting_state/vein/native/vein.nii",
        ]

input_surf = [
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/lh.layer0",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/lh.layer1",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/lh.layer2",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/lh.layer3",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/lh.layer4",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/lh.layer5",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/lh.layer6",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/lh.layer7",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/lh.layer8",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/lh.layer9",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/lh.layer10",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/rh.layer0",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/rh.layer1",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/rh.layer2",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/rh.layer3",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/rh.layer4",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/rh.layer5",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/rh.layer6",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/rh.layer7",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/rh.layer8",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/rh.layer9",
        "/data/pt_01880/Experiment2_Rivalry/p3/anatomy/layer/rh.layer10",
        ]

# parameters
deformation = "/data/pt_01880/Experiment2_Rivalry/p3/deformation/resting_state/epi2orig.nii.gz"
path_output = "/data/pt_01880/Experiment2_Rivalry/p3/resting_state/vein"
interpolation = 'linear' # can be linear or nearest

""" do not edit below """

# apply deformation    
path_def = os.path.join(path_output,"def")
for i in range(len(input_file)):
    apply_coordinate_mappings(input_file[i], 
                              deformation, 
                              interpolation=interpolation, 
                              padding='closest', 
                              save_data=True, 
                              overwrite=True, 
                              output_dir=path_def,
                              file_name=None,
                              )

# map to ana
path_surf = os.path.join(path_output,"surf")
for i in range(len(input_file)):
    # name of input volume
    filename_def = os.path.splitext((os.path.basename(input_file[i])))[0]

    for j in range(len(input_surf)):
        # hemisphere
        hemi = os.path.splitext(os.path.basename(input_surf[j]))[0]
    
        # sample on surface
        map2surface(input_surf[j], 
                    os.path.join(path_def,filename_def+"_def-img.nii.gz"), 
                    hemi, 
                    path_surf, 
                    input_white=None, 
                    input_ind=None, 
                    cleanup=True)
