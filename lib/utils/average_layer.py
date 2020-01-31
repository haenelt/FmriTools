def average_layer(img_input, path_output, basename_output, mode="mean"):
    """
    This averages data across different layers or sessions. Input arrays should be in mgh format. 
    The output gets the suffix of the chosen mode.
    Inputs:
        *img_input: list of input layers.
        *path_output: path where output is saved.
        *basename_output: basename of written output file.
        *mode: average mode (mean or median).
    
    created by Daniel Haenelt
    Date created: 25-10-2019
    Last modified: 31-01-2020
    """
    import sys
    import os
    import numpy as np
    import nibabel as nb
    
    # initialise array
    data = nb.load(img_input[0])
    data_size = data.header["dims"][0]
    data_res = np.zeros((data_size,len(img_input)))
    
    # collect input arrays
    for i in range(len(img_input)):
        data_res[:,i] = nb.load(img_input[i]).get_fdata()[:,0,0]
    
    # average
    if mode is "mean":
        data_res = np.mean(data_res, axis=1)
    elif mode is "median":
        data_res = np.median(data_res, axis=1)
    else:
        sys.exit("Choose a valid mode!")
    
    # expand dimensions to match with input array
    data_res = np.expand_dims(data_res, axis=1)
    data_res = np.expand_dims(data_res, axis=1)
    
    # write output file
    output = nb.Nifti1Image(data_res, data.affine, data.header)
    nb.save(output,os.path.join(path_output, basename_output+"_"+mode+".mgh"))