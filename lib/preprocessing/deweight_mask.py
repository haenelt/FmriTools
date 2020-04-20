def deweight_mask(file_in, mask_in, mask_max=0.25, sigma_gaussian=10.0, write_output=None, 
                  path_output=None):
    """
    This function computes a binary mask by pooling all voxels above a given threshold and replaces
    all image voxels by its gaussian filtered image voxels within this binary mask.
    Inputs:
        *file_in: filename of input image.
        *mask_in: filename of input mask.
        *mask_max: cutoff threshold.
        *sigma_gaussian: sigma for gaussian filter.
        *write_output: write output image (boolean).
        *path_output: path where output is written.
    Outputs:
        *data_array: image matrix with filtered voxels.
        
    created by Daniel Haenelt
    Date created: 20-04-2020 
    Last modified: 20-04-2020
    """
    import os
    import nibabel as nb
    from scipy.ndimage import gaussian_filter
    from lib.io.get_filename import get_filename
    
    # get basename of phase file
    _, name_file, ext_file = get_filename(file_in)
    
    # load unwrapped phase data
    data = nb.load(file_in)
    data_array = data.get_fdata()
    
    # load standard deviation data
    mask = nb.load(mask_in)
    mask_array = mask.get_fdata()

    # threshold standard deviation    
    mask_array[mask_array < mask_max] = 0
    mask_array[mask_array != 0] = 1

    # apply gaussian filter to phase data
    data_array_gaussian = gaussian_filter(data_array, 
                                          sigma_gaussian, 
                                          order=0, 
                                          output=None, 
                                          mode='reflect', 
                                          cval=0.0, 
                                          truncate=4.0)
    
    # replace data
    data_array[mask_array == 1] = data_array_gaussian[mask_array == 1]
    
    # write output
    if write_output:
        output = nb.Nifti1Image(data_array, data.affine, data.header)
        nb.save(output,os.path.join(path_output,name_file+"_filtered"+ext_file))
    
    return data_array