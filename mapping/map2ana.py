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
file_in = ["/nobackup/actinium2/haenelt/V2STRIPES/p2/colour/results/spmT/def/spmT_colour_bw_GE_EPI1_def.nii.gz",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/colour/results/spmT/def/spmT_colour_bw_GE_EPI2_def.nii.gz",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/colour/results/spmT/def/spmT_colour_bw_GE_EPI_def_mean.nii.gz",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/colour/results/spmT/def/spmT_colour_bw_SE_EPI1_def.nii.gz",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/colour/results/spmT/def/spmT_colour_bw_SE_EPI2_def.nii.gz",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/colour/results/spmT/def/spmT_colour_bw_SE_EPI_def_mean.nii.gz",
           ]

surf_in = ["/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer0",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer1",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer2",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer3",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer4",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer5",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer6",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer7",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer8",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer9",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer10",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer11",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer12",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer13",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer14",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer15",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer16",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer17",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer18",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer19",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/lh.layer20",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer0",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer1",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer2",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer3",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer4",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer5",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer6",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer7",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer8",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer9",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer10",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer11",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer12",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer13",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer14",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer15",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer16",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer17",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer18",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer19",
           "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab/layer/rh.layer20",
           ]

# parameters
deformation_in1 = None
deformation_in2 = None
path_output = "/nobackup/actinium2/haenelt/V2STRIPES/p2/sab"
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