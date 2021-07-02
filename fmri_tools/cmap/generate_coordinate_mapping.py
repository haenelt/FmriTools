# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb
from numpy.matlib import repmat


def generate_coordinate_mapping(file_in, pad, path_output=None, suffix=None,
                                time=False, write_output=False):
    """Generate coordinate mapping.
    
    Generates coordinate mapping for an input volume. Either one or multiple 
    coordinate maps are saved in the output folder depending on the 
    dimensionality (3d or 4d) of the input image. Image padding can be applied 
    which expands the image matrix of each axis in both directions.    

    Parameters
    ----------
    file_in : str
        File name of input file.
    pad : int
        Image padding size.
    path_output : str, optional
        Path where output is saved. The default is None.
    suffix : str, optional
        Add suffix to file name. The default is None.
    time : bool, optional
        Compute coordinate map for each time step. The default is False.
    write_output : bool, optional
        Write nifti volume. The default is False.

    Returns
    -------
    output : niimg
        Coordinate mapping.
    
    """

    # create output folder
    if path_output:
        if not os.path.exists(path_output):
            os.makedirs(path_output)

    # load data
    data_img = nb.load(file_in)

    # get matrix size
    x_size = data_img.header["dim"][1] + 2 * pad
    y_size = data_img.header["dim"][2] + 2 * pad
    z_size = data_img.header["dim"][3] + 2 * pad

    if time is False:
        t_size = 1
    else:
        t_size = data_img.header["dim"][4]

    # define coordinate
    coordinate_mapping = np.zeros((x_size, y_size, z_size, 3), dtype='float')

    # coordinate mapping in x-direction
    x = np.array(np.arange(-pad, x_size - pad, 1), dtype='float')
    x = np.transpose(repmat(x, y_size, 1))
    x = np.dstack([x] * z_size)

    # coordinate mapping in y-direction
    y = np.array(np.arange(-pad, y_size - pad), dtype='float')
    y = repmat(y, x_size, 1)
    y = np.dstack([y] * z_size)

    # coordinate mapping in z-direction
    z = np.ones((x_size, y_size, z_size))
    z = np.arange(-pad, z_size - pad) * z

    # merge directions
    coordinate_mapping[:, :, :, 0] = x
    coordinate_mapping[:, :, :, 1] = y
    coordinate_mapping[:, :, :, 2] = z

    # write coordinate mapping
    output = nb.Nifti1Image(coordinate_mapping, data_img.affine,
                            nb.Nifti1Header())

    # write coordinate mapping for each time point   
    if t_size == 1:
        if write_output:
            file_out = os.path.join(path_output, 'cmap_' + suffix + '.nii')
            nb.save(output, file_out)
    else:
        for i in range(t_size):
            if write_output:
                file_out = os.path.join(path_output,
                                        'cmap_' + suffix + '_' + str(i) +
                                        '.nii')
                nb.save(output, file_out)

    return output
