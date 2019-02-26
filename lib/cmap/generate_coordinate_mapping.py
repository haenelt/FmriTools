def generate_coordinate_mapping(input, pad, path_output, suffix="", time=False):
    """
    Generates coordinate mapping for an input volume. Either one or multiple coordinate maps are 
    saved in the output folder depending on the dimensionality (3d or 4d) of the input image. Image
    padding can be applied which expands the image matrix of each axis in both directions.
    Inputs:
        *input: input file.
        *pad: image padding size.
        *path_output: path where output is saved.
        *suffix: add suffix to file name.
        *time: compute coordinate map for each time step.

    created by Daniel Haenelt
    Date created: 21-11-2018             
    Last modified: 31-01-2019
    """
    import os
    import numpy as np
    import nibabel as nb
    from numpy.matlib import repmat
    
    # create output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    
    # load data
    data_img = nb.load(input)
    
    # get matrix size
    x_size = data_img.header["dim"][1] + 2*pad
    y_size = data_img.header["dim"][2] + 2*pad
    z_size = data_img.header["dim"][3] + 2*pad
    
    if time is False:
        t_size = 1
    else:
        t_size = data_img.header["dim"][4]

    # define coordinate
    coordinate_mapping = np.zeros((x_size, y_size, z_size, 3), dtype='float')

    # coordinate mapping in x-direction
    X = np.array(np.arange(-pad,x_size-pad,1), dtype='float')
    X = np.transpose(repmat(X, y_size, 1))
    X = np.dstack([X]*z_size)
    
    # coordinate mapping in y-direction
    Y = np.array(np.arange(-pad,y_size-pad), dtype='float')
    Y = repmat(Y, x_size, 1)
    Y = np.dstack([Y]*z_size)
    
    # coordinate mapping in z-direction
    Z = np.ones((x_size, y_size, z_size))
    Z = np.arange(-pad,z_size-pad) * Z
    
    # merge directions
    coordinate_mapping[:,:,:,0] = X
    coordinate_mapping[:,:,:,1] = Y
    coordinate_mapping[:,:,:,2] = Z

    # write coordinate mapping
    output = nb.Nifti1Image(coordinate_mapping, data_img.affine, data_img.header)
    output.set_data_dtype(np.float)

    # write coordinate mapping for each time point   
    if t_size == 1:
        fileOUT = os.path.join(path_output,'cmap_'+suffix+'.nii')
        nb.save(output,fileOUT)
    else:
        for i in range(t_size):
            fileOUT = os.path.join(path_output,'cmap_'+suffix+'_'+str(i)+'.nii')
            nb.save(output,fileOUT)

    