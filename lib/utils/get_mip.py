def get_mip(input_file, slice_start, slice_end, path_output):
    """
    Create minimum intensity projection between between slices of a BOLD time series.
    Inputs:
        *input_file: input image.
        *slice_start: starting slice index for minimum intensity projection.
        *slice_end: ending slice index for minimum intensity projection.
        *path_output: path where to write output image.
        
    created by Daniel Haenelt
    Date created: 16-11-2018
    Last modified: 16-11-2018
    """
    import os
    import numpy as np
    import nibabel as nb
    
    # make subfolders
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    
    # load data
    data_img = nb.load(input_file)
    data_array = data_img.get_data()
    
    # mean epi
    mean_array = np.mean(data_array, axis=3)
    
    # new image size
    nx = data_img.shape[1]
    ny = data_img.shape[2]
    
    # get minimum projection
    slice_array = mean_array[:, :, slice_start:slice_end]
    
    # mip
    mip = np.zeros([nx, ny])
    for i in range(nx):
        for j in range(ny):
            mip[i, j] = np.min(slice_array[i, j])
    
    # output image matrix
    newimg = nb.Nifti1Image(mip, data_img.affine, data_img.header)
    newimg.header["dim"][0] = 2
    newimg.header["dim"][1] = nx
    newimg.header["dim"][2] = ny
    newimg.header["dim"][3] = 1
    newimg.header["dim"][4] = 1
    nb.save(newimg, os.path.join(path_output,"mip.nii"))
