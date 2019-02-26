def get_gaussian(file, path_output, mu=0, sigma=1):
    """
    Creates a NIfTI volume with Gaussian noise as voxel intensities and saves a nifti volume with 
    matrix size and header information taken from the input image.
    Inputs:
        *file: input image.
        *mu: mean of Gaussian noise.
        *sigma: standard deviation of Gaussian noise.
        *path_output: path where output file is written.
    
    created by Daniel Haenelt
    Date created: 20-11-2018
    Last modified: 20-11-2018
    """
    import os
    import numpy as np
    import nibabel as nb
    
    # load img
    img = nb.load(input_file)
    
    # new image size
    x_size = img.shape[0]
    y_size = img.shape[1]
    z_size = img.shape[2]
    
    # random Gaussian noise for each pixel
    img_data = np.random.normal(mu,sigma,x_size*y_size*z_size)
    img_data = np.reshape(img_data,(x_size,y_size,z_size))
    
    # write output image
    newimg = nb.Nifti1Image(img_data, img.affine, img.header)
    nb.save(newimg,os.path.join(path_output, "gaussian_noise.nii"))
