def map2grid(file_grid, file_input, sigma, path_output, overwrite=True):
    """
    This script allows you to sample indexed morphological data onto the regular grid. Optional, a 
    gaussian filter can be applied to the output image.
    Inputs:
        *file_grid: filename of grid coordinate mapping.
        *file_input: filename of morphological data or *.mgh data.
        *sigma: standard deviation of Gaussian kernel.
        *path_output: path where output is saved.
    Output:
        *grid_array: file mapped onto array.
    
    created by Daniel Haenelt
    Date created: 01-11-2018             
    Last modified: 15-02-2019
    """
    import os
    import numpy as np
    import nibabel as nb
    from nibabel.freesurfer.io import read_morph_data
    from scipy.ndimage.filters import gaussian_filter

    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # load data
    grid_img = nb.load(file_grid)
    grid_array = grid_img.get_fdata()
    if os.path.splitext(file_input)[1] == ".mgh":
        morph = nb.load(file_input).get_fdata()
    else:
        morph = read_morph_data(file_input)

    # sample data onto grid
    for i in range(np.size(grid_array,0)):
        for j in range(np.size(grid_array,1)):
            if grid_array[i,j] != 0:
                grid_array[i,j] = morph[grid_array[i,j].astype(int)]
    
    # gaussian filter (opt)
    if sigma != 0:
        order = 0
        mode = "reflect"
        truncate = 4.0
        grid_array = gaussian_filter(grid_array,
                                     sigma=sigma,
                                     order=order,
                                     mode=mode,
                                     truncate=truncate)
    
    # write output data
    if overwrite:
        if sigma != 0:
            filenameOUT = os.path.join(path_output,os.path.splitext(os.path.basename(file_input))[0]+"_grid.nii")
        else:
            filenameOUT = os.path.join(path_output,os.path.splitext(os.path.basename(file_input)+"_sigma"+str(sigma))[0]+"_grid.nii")
        
        output = nb.Nifti1Image(grid_array, grid_img.affine, grid_img.header)
        nb.save(output,filenameOUT)

    return grid_array