"""
Map2ana

The purpose of the following script is to sample a registered volume to its corresponding surface
mesh.

created by Daniel Haenelt
Date created: 06-03-2019            
Last modified: 06-03-2019  
"""

import os
from nighres.registration import apply_coordinate_mappings
from lib.mapping import map2surface

# input
input_file = ["/data/pt_01880/V2STRIPES/p2/colour/results/spmT/native/spmT_colour_bw_GE_EPI2.nii",
              ]

input_surf = ["/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/lh.layer0",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/lh.layer1",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/lh.layer2",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/lh.layer3",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/lh.layer4",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/lh.layer5",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/lh.layer6",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/lh.layer7",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/lh.layer8",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/lh.layer9",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/rh.layer0",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/rh.layer1",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/rh.layer2",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/rh.layer3",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/rh.layer4",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/rh.layer5",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/rh.layer6",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/rh.layer7",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/rh.layer8",
	      "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/layer/rh.layer9",
        ]

deformation = "/data/pt_01880/V2STRIPES/p2/deformation/GE_EPI2/epi2ana.nii.gz"
path_output = "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p2/color/ge_epi2"

""" do not edit below """

# apply deformation
path_def = os.path.join(path_output,"def")
for i in range(len(input_file)):
    apply_coordinate_mappings(input_file[i], 
                              deformation, 
                              interpolation='linear', 
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
