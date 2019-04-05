def get_tsnr(input, path_output, name_output, TR, cutoff_highpass, pathSPM):
    """
    This function computes the tsnr of one time series.
    Inputs:
        *input: input time series.
        *path_output: path where to save mean image
        *name_output: output file name without file extension.
        *TR: repetition time in s.
        *cutoff_highpass: cutoff in s for baseline correction.
        *pathSPM: path to SPM12.
        
    created by Daniel Haenelt
    Date created: 05-02-2019         
    Last modified: 05-02-2019
    """
    import os
    import numpy as np
    import nibabel as nb

    # make subfolders
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # change to lib folder
    path_baseline_func = os.path.realpath(__file__)
    path_baseline_func = os.path.dirname(os.path.dirname(path_baseline_func))
    path_baseline_func = os.path.join(path_baseline_func,"preprocessing")
    os.chdir(path_baseline_func)
    
    # prepare path and filename
    path = os.path.split(input)[0]
    file = os.path.split(input)[1]

    # look for baseline corrected time series
    if not os.path.isfile(os.path.join(path,"b"+file)):
        os.system("matlab" + \
                  " -nodisplay -nodesktop -r " + \
                  "\"baseline_correction(\'{0}\', {1}, {2}, \'{3}\'); exit;\"". \
                  format(input, TR, cutoff_highpass, pathSPM))

    # load baseline corrected time series
    data_img = nb.load(os.path.join(path,"b"+file))
    data_array = data_img.get_fdata()
    
    # get mean and std
    data_mean_array = np.mean(data_array, axis=3)
    data_std_array = np.std(data_array, axis=3)
    data_std_array[data_std_array == 0] = np.nan # set zeroes to nan
    
    # get tsnr of time series
    data_tsnr_array = data_mean_array / data_std_array
    data_tsnr_array[np.isnan(data_tsnr_array)] = 0
    
    # write output
    data_img.header["dim"][0] = 3
    data_img.header["dim"][4] = 1
    data_img.header["pixdim"][3] = 1
    
    data_img = nb.Nifti1Image(data_tsnr_array, data_img.affine, data_img.header)
    nb.save(data_img, os.path.join(path_output,"tsnr"+name_output+".nii"))
