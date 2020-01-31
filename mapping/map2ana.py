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
file_in = [        
        "/data/pt_01880/Experiment3_Stripes/p3/mpm/pd_kp_mtflash3d_v1ax_0p5_0008/Results/s2803151-124913-00001-00352-1_PD_gnlcorr_scaled_inverted.nii",
        "/data/pt_01880/Experiment3_Stripes/p3/mpm/pd_kp_mtflash3d_v1ax_0p5_0008/Results/s2803151-124913-00001-00352-1_R1_gnlcorr.nii",
        "/data/pt_01880/Experiment3_Stripes/p3/mpm/pd_kp_mtflash3d_v1ax_0p5_0008/Results/s2803151-124913-00001-00352-1_R2s_WOLS_gnlcorr_scaled.nii",
        ]

surf_in = [
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer0_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer1_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer2_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer3_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer4_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer5_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer6_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer7_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer8_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer9_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer10_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer11_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer12_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer13_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer14_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer15_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer16_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer17_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer18_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer19_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/lh.layer20_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer0_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer1_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer2_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer3_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer4_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer5_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer6_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer7_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer8_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer9_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer10_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer11_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer12_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer13_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer14_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer15_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer16_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer17_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer18_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer19_mpm",
        "/data/pt_01880/Experiment3_Stripes/p3/anatomy/layer_mpm/rh.layer20_mpm",
        ]

# parameters
deformation_in = None
path_output = "/data/pt_01880/Experiment3_Stripes/p3/mpm/pd_kp_mtflash3d_v1ax_0p5_0008/Results"
interpolation = "linear" # can be linear or nearest

""" do not edit below """

# output folders
path_def = os.path.join(path_output,"def")
path_surf = os.path.join(path_output,"surf")

for i in range(len(file_in)):
    
    # apply deformation    
    if deformation_in is not None:
        apply_coordinate_mappings(file_in[i], 
                                  deformation_in, 
                                  interpolation=interpolation, 
                                  padding='closest', 
                                  save_data=True, 
                                  overwrite=True, 
                                  output_dir=path_def,
                                  file_name=None,
                                  )
    
    # get filename of deformed file
    if deformation_in is not None:
        _, name_file, _ = get_filename(file_in[i])
        filename_def = os.path.join(path_def, name_file+"_def-img.nii.gz") 
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