# -*- coding: utf-8 -*-

# python standard library inputs
import sys

# external inputs
import numpy as np


def read_patch(file_in):
    """Read patch.

    This function reads an freesurfer patch saved in binary format. It is an 
    equivalent to the matlab function read_patch in the ./freesurfer/matlab 
    folder. Data is read in big endian order.    

    Parameters
    ----------
    file_in : str
        Full path of the input file.

    Returns
    -------
    x : ndarray
        x-coordinates of patch.
    y : ndarray
        y-coordinates of patch.
    z : ndarray
        z-coordinates of patch (if flattened, there should be only zeros).
    ind : ndarray
        index of each node.
    
    """

    # load data
    data_array_int = np.fromfile(file_in, np.dtype(">i"))
    data_array_float = np.fromfile(file_in, np.dtype(">f"))

    # check version
    ver = data_array_int[0]
    if not ver == -1:
        sys.exit("incorrect version # " + str(ver) + " (not -1) found in file")

    # size of data array
    data_size = data_array_int[1]

    # reshape data into array
    data_array_int = data_array_int[2:]
    data_array_float = data_array_float[2:]
    data_array_int = data_array_int.reshape(data_size, 4)
    data_array_float = data_array_float.reshape(data_size, 4)

    # get vertex indices and coordinates
    ind = data_array_int[:, 0]
    ind[ind < 0] = -ind[ind < 0] - 1
    ind[ind >= 0] = ind[ind >= 0] - 1

    x = data_array_float[:, 1]
    y = data_array_float[:, 2]
    z = data_array_float[:, 3]

    return x, y, z, ind
