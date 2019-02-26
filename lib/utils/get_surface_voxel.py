def get_surface_voxel(input_file, path_output):
    """
    Compute surface voxels for to illustrate curvature dependence in surface representations as 
    introduced in Kay et al. (2018).
    Inputs:
        *input_file: input image.
        *path_output: path where to write output image.
    
    created by Daniel Haenelt
    Date created: 19-11-2018
    Last modified: 19-11-2018
    """
    import os
    import numpy as np
    import nibabel as nb

    # parameters
    x_calc = 1
    y_calc = 1
    z_calc = 1
    x_max = 1
    y_max = 2
    z_max = 4

    # make subfolders
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    
    # load img
    img = nb.load(input_file)
    
    # img range
    x_size = img.shape[0]
    y_size = img.shape[1]
    z_size = img.shape[2]
    
    # compute surface voxels
    img_res = np.zeros([x_size,y_size,z_size])
    for i in range(x_size):
        for j in range(y_size):
            for k in range(z_size):
                if np.mod(i,2) != 0 and x_calc != 0:
                    img_res[i,j,k] += x_max
                if np.mod(j,2) != 0 and y_calc != 0:
                    img_res[i,j,k] += y_max
                if np.mod(k,2) != 0 and z_calc != 0:
                    img_res[i,j,k] += z_max    
             
    # write output image
    newimg = nb.Nifti1Image(img_res, img.affine, img.header)
    nb.save(newimg,os.path.join(path_output,"surface_voxel.nii"))
