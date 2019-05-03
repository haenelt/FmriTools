def average_time_series(img_input, path_output, name_output):
    """
    This function computes the element-wise average of multiple nifti files.
    Inputs:
        *img_input: list of nifti input paths.
        *path_output: path where output is saved.
        *name_output: basename of output files.
    
    created by Daniel Haenelt
    Date created: 03-05-2019
    Last modified: 03-05-2019
    """
    import os
    import numpy as np
    import nibabel as nb

    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # load first dataset to initialize final time series
    res = nb.load(img_input[0])
    res_array = np.zeros_like(res.get_fdata())
    
    for i in range(len(img_input)):    
        res_array += nb.load(img_input[i]).get_fdata()
    
    # divide summed time series by number of time series
    res_array = res_array / len(img_input)
    
    # write output
    output = nb.Nifti1Image(res_array, res.affine, res.header)
    nb.save(output,os.path.join(path_output,name_output+".nii"))
    