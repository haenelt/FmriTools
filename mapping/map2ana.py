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
input_file = ["/nobackup/eminem2/attar/3T_Connectom/pilot/subj-03/fmri/retinotopy/avg/native/ecc_imag_avg.nii",
              "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-03/fmri/retinotopy/avg/native/ecc_phase_avg.nii",
              "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-03/fmri/retinotopy/avg/native/ecc_real_avg.nii",
              "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-03/fmri/retinotopy/avg/native/pol_imag_avg.nii",
              "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-03/fmri/retinotopy/avg/native/pol_phase_avg.nii",
              "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-03/fmri/retinotopy/avg/native/pol_real_avg.nii",
              ]

input_surf = ["/nobackup/eminem2/attar/3T_Connectom/pilot/subj-03/smri/mp2rage/2.0-analysis/5.0-surface/freesurfer/surf/lh.mid",
	      "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-03/smri/mp2rage/2.0-analysis/5.0-surface/freesurfer/surf/rh.mid",
        ]

deformation = "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-03/fmri/retinotopy/deformation/epi2orig.nii.gz"
path_output = "/nobackup/eminem2/attar/3T_Connectom/pilot/subj-03/fmri/retinotopy/avg"

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
