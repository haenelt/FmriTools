"""
Transform time series to target space

In the following script, epi time series in native space are transformed to a target space using a
deformation field. The transformed time series get the prefix r.

created by Daniel Haenelt
Date created: 07-08-2019            
Last modified: 07-08-2019  
"""
import os
from sh import gunzip
from nighres.registration import apply_coordinate_mappings

# input
input_epi = [
        "/home/raid2/haenelt/Desktop/udata_test.nii",
        ]

input_reg = [
        "/home/raid2/haenelt/Desktop/epi2orig.nii.gz",
        ]

# parameters
interpolation = "linear"
padding = "closest"

""" do not edit below """

# apply deformation    
if len(input_epi) == len(input_reg):
    for i in range(len(input_epi)):
        apply_coordinate_mappings(input_epi[i], 
                                  input_reg[i], 
                                  interpolation=interpolation, 
                                  padding=padding,
                                  save_data=True, 
                                  overwrite=True, 
                                  output_dir=os.path.dirname(input_epi[i]),
                                  file_name=None,
                                  )

        # time series path and basename
        path = os.path.dirname(input_epi[i])
        file = os.path.splitext(os.path.basename(input_epi[i]))[0]

        # unzip output
        gunzip(os.path.join(path,file+"_def-img.nii.gz"))
        
        # rename output time series
        os.rename(os.path.join(path,file+"_def-img.nii"),os.path.join(path,"r"+file+".nii"))

else:
    print("Number of time series and deformation are not the same!")