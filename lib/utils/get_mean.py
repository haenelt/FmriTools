def get_mean(input, path_output, name_output, type="mean"):
    """
    This function computes the mean image of one or more time series.
    Inputs:
        *input: single file or list of files.
        *path_output: path where to save mean image
        *name_output: output file name without file extension.
        *type: can be either mean or median.
        
    created by Daniel Haenelt
    Date created: 04-02-2019         
    Last modified: 24-10-2019
    """
    import os
    import copy
    import numpy as np
    import nibabel as nb
    
    # make subfolders
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    
    if len(np.shape(input)) > 0:
        
        # get dimensions
        data_img = nb.load(input[0])
        
        x_size = data_img.header["dim"][1]
        y_size = data_img.header["dim"][2]
        z_size = data_img.header["dim"][3]
        
        t_size = 0
        for i in range(len(input)):
            data_img = nb.load(input[i])
            t_size += data_img.header["dim"][4]
        
        data_all_array = np.zeros([x_size, y_size, z_size, t_size])
        t_start = 0
        for i in range(len(input)):
            data_img = nb.load(input[i])
            data_array = data_img.get_fdata()
            t_end = data_img.header["dim"][4] + t_start     
            data_all_array[:,:,:,t_start:t_end] = data_array
            t_start = copy.deepcopy(t_end)
                
    else:
               
        # load data
        data_img = nb.load(input)
        data_all_array = data_img.get_fdata()
    
    # calculate mean
    if type is "mean":
        data_mean_array = np.mean(data_all_array,axis=3)
    elif type is "median":
        data_mean_array = np.median(data_all_array,axis=3)
    else:
        print("Choose a valid mean type!")
    
    # write output
    if len(np.shape(input)) > 0:
        data_img = nb.load(input[0])
    else:
        data_img = nb.load(input)
    data_img.header["dim"][0] = 3
    data_img.header["dim"][4] = 1
    data_img.header["pixdim"][3] = 1
    
    # write mean image
    mean_img = nb.Nifti1Image(data_mean_array, data_img.affine, data_img.header)
    nb.save(mean_img,os.path.join(path_output,"mean_"+name_output+".nii"))