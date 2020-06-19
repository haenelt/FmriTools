def remove_nans(file_in, file_out):
    """
    Reads a nifti volume and sets all nans to zero. 
    Inputs:
        *file_in: filename of input volume.
        *file_out: filename of output volume.
        
    created by Daniel Haenelt
    Date created: 19-06-2020             
    Last modified: 19-06-2020
    """
    import numpy as np
    import nibabel as nb

    # load data
    data_img = nb.load(file_in)
    data_array = data_img.get_fdata()
    
    # set nans to zero
    data_array[np.isnan(data_array)] = 0
    
    # write output data
    output = nb.Nifti1Image(data_array, data_img.affine, data_img.header)
    nb.save(output, file_out)