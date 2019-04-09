def get_laminar_profile(input, path_output, hemi, name_output, mode):
    """
    This function computes the laminar profile (mean, median, max, min) of data sampled onto
    different cortical layers.
    Inputs:
        *input: List of sampled data onto different cortical layers (surface mesh).
        *path_output: path where to save laminar profile.
        *hemi: hemisphere.
        *name_output: output file name without file extension.
        *mode: mean, median, max, min.
        
    created by Daniel Haenelt
    Date created: 09-04-2019         
    Last modified: 09-04-2019
    """
    import os
    import numpy as np
    import nibabel as nb

    # make subfolders
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # get input header
    data_img = nb.load(input[0])

    # read input
    data = []
    for i in range(len(input)):
        data.append(np.squeeze(nb.load(input[i]).get_fdata()))
        
    # convert input into array
    data = np.array(data)

    # compute laminar profile
    if mode == "mean":
        data = np.mean(data,0)
    elif mode == "median":
        data = np.median(data,0)
    elif mode == "max":
        data = np.max(data,0)
    elif mode == "min":
        data = np.min(data,0)
    else:
        print("Choose a valid mode!")

    # expand axes to match input array
    data = np.expand_dims(data,1)
    data = np.expand_dims(data,1)

    # write output image
    output = nb.Nifti1Image(data, data_img.affine, data_img.header)
    nb.save(output, os.path.join(path_output,hemi+"."+name_output+".mgh"))
