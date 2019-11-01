def get_mean4d(input, path_output, name_output):
    """
    This function computes the mean time series of one or more time series.
    Inputs:
        *input: list of 4d nifti files.
        *path_output: path where to save mean image
        *name_output: output file name without file extension.
        
    created by Daniel Haenelt
    Date created: 31-10-2019         
    Last modified: 31-10-2019
    """
    import os
    import numpy as np
    import nibabel as nb
    
    # make subfolders
    if not os.path.exists(path_output):
        os.makedirs(path_output)
        
    # get dimensions
    data_img = nb.load(input[0])
    
    res_array = np.zeros_like(data_img.get_fdata())
    for i in range(len(input)):
        res_array += nb.load(input[i]).get_fdata()

    res_array = res_array / len(input)
    
    # write mean time series
    mean_img = nb.Nifti1Image(res_array, data_img.affine, data_img.header)
    nb.save(mean_img,os.path.join(path_output,"mean_"+name_output+".nii"))