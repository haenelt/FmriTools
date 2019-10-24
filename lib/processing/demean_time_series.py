def demean_time_series(img_input):
    """
    This function demeans each voxel time series. Input is either a 4d nifti or compressed nifti 
    file.
    Inputs:
        *img_input: 4d nifti volume.
    
    created by Daniel Haenelt
    Date created: 24-10-2019
    Last modified: 24-10-2019
    """
    import os
    import numpy as np
    import nibabel as nb

    # get path and filename for output
    path = os.path.dirname(img_input)
    if os.path.splitext(img_input)[1] == ".gz":
        file = os.path.splitext(os.path.splitext(os.path.basename(img_input))[0])[0]
    else:
        file = os.path.splitext(os.path.basename(img_input))[0]

    # load first dataset to initialize final time series
    data = nb.load(img_input)
    data_array = data.get_fdata()
    
    # get mean of each voxel time series
    data_mean = np.mean(data_array, axis=3)
    
    # demean time series
    for i in range(np.shape(data_array)[3]):
        data_array[:,:,:,i] = ( data_array[:,:,:,i] - data_mean )/ data_mean * 100
       
    # write output
    output = nb.Nifti1Image(data_array, data.affine, data.header)
    if os.path.splitext(img_input)[1] == ".gz":
        nb.save(output,os.path.join(path,file+"_demean.nii.gz"))
    else:
        nb.save(output,os.path.join(path,file+"_demean.nii"))
    