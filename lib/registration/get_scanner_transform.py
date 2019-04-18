def get_scanner_transform(input_source, input_target, path_output):
    """
    This function uses the scanner coordinates to create a coordinate map between two images in the 
    same scanner coordinate system. The orientation matrices written in the header of both files 
    are taken to get the trasformation between both images. The output contains a 4d coordinate map 
    describing the transformation from source to target image. Input files should be in nifti 
    format.
    Inputs:
        *input_source: absolute path of source file.
        *target_source: absolute path of target file.
        *path_output: path where output is saved.

    created by Daniel Haenelt
    Date created: 13-11-2018       
    Last modified: 22-01-2019
    """
    import os
    import numpy as np
    import nibabel as nb
    import numpy.linalg as npl
    from nibabel.affines import apply_affine
    from numpy.matlib import repmat
       
    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)
    
    # load source and target images
    source_img = nb.load(input_source)
    target_img = nb.load(input_target)
    
    # get affine transformation
    target2source = npl.inv(source_img.affine).dot(target_img.affine)
    
    # initialise target coordinate map
    x_size = target_img.header["dim"][1]
    y_size = target_img.header["dim"][2]
    z_size = target_img.header["dim"][3]
    
    coordinate_mapping = np.zeros((x_size,y_size,z_size,3), dtype='float')
    
    # coordinate mapping in x-direction
    X = np.array(np.arange(0,x_size,1), dtype='float')
    X = np.transpose(repmat(X, y_size, 1))
    X = np.dstack([X]*z_size)
    
    # coordinate mapping in y-direction
    Y = np.array(np.arange(0,y_size), dtype='float')
    Y = repmat(Y, x_size, 1)
    Y = np.dstack([Y]*z_size)
    
    # coordinate mapping in z-direction
    Z = np.ones((x_size, y_size, z_size))
    Z = np.arange(0,z_size) * Z
    
    # merge directions
    coordinate_mapping[:,:,:,0] = X
    coordinate_mapping[:,:,:,1] = Y
    coordinate_mapping[:,:,:,2] = Z
    
    # apply transformation to coordinates
    coordinate_mapping = apply_affine(target2source,coordinate_mapping)
    
    # write coordinate map
    target_img.header["dim"][0] = 4
    target_img.header["dim"][4] = 3
    target_img.set_data_dtype(np.float)
       
    # get filenames
    if os.path.splitext(os.path.basename(input_source))[1] == '.gz':
        name_source = os.path.splitext(os.path.splitext(os.path.basename(input_source))[0])[0]
    else:
        name_source = os.path.splitext(os.path.basename(input_source))[0]
    
    if os.path.splitext(os.path.basename(input_target))[1] == '.gz':
        name_target = os.path.splitext(os.path.splitext(os.path.basename(input_target))[0])[0]
    else:
        name_target = os.path.splitext(os.path.basename(input_target))[0]
      
    output = nb.Nifti1Image(coordinate_mapping, target_img.affine, target_img.header)
    nb.save(output,os.path.join(path_output,name_source+"_2_"+name_target+"_scanner.nii"))