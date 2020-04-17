"""
MPM scaling

This scripts restricts values in a nifti volume to a defined data range. Ceiling outside are set to
the set maximum or minimum value. Optionally, the range is inverted.

created by Daniel Haenelt
Date created: 03-03-2020
Last modified: 03-03-2020
"""
import os
import numpy as np
import nibabel as nb
from lib.io.get_filename import get_filename

# input
input = [
    "/data/pt_01983/anat/Session1/mpms/afi/Results/s2019-10-17_13-17-140642-00001-00352-1__dis3d_R1.nii",
    "/data/pt_01983/anat/Session1/mpms/seste/Results/s2019-10-17_13-17-140642-00001-00352-1__dis3d_R1.nii",
    "/data/pt_01983/anat/Session2/mpms/afi/Results/s2019-10-22_12-57-140400-00001-00352-1__dis3d_R1.nii",
    "/data/pt_01983/anat/Session2/mpms/seste/Results/s2019-10-22_12-57-140400-00001-00352-1__dis3d_R1.nii",
    ]

val_min = 0 # if not None
val_max = 2 # if not None
val_invert = False # if True

""" do not edit below """

for i in range(len(input)):
    
    # get filename
    path_file, name_file, ext_file = get_filename(input[i])
    
    # load data
    data = nb.load(input[i])
    data_array = data.get_fdata()
    
    # scale
    if val_min is not None:
        data_array[data_array < val_min] = val_min
        
    if val_max is not None:
        data_array[data_array > val_max] = val_max

    # invert    
    if val_invert:
        data_array = np.max(data_array) - data_array
    
    # write output
    output = nb.Nifti1Image(data_array, data.affine, data.header)
    if val_invert:
        file_out =  os.path.join(path_file, name_file + "_scaled_inverted" + ext_file)
    else:
        file_out =  os.path.join(path_file, name_file + "_scaled" + ext_file)
    
    nb.save(output, file_out)