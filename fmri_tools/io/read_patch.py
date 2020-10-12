# -*- coding: utf-8 -*-

# python standard library inputs
import sys

# external inputs
import numpy as np


def read_patch(filename):
    """
    This function reads an freesurfer patch saved in binary format. It is an 
    equivalent to the matlab function read_patch in the ./freesurfer/matlab 
    folder. Data is read in big endian order.
    Inputs:
        *filename: full path of the input file.
    Outputs:
        *x: x-coordinates of patch.
        *y: y-coordinates of patch.
        *z: z-coordinates of patch (if flattened, there should be only zeros).
        *ind: index of easch node.
    
    created by Daniel Haenelt
    Date created: 01-11-2018             
    Last modified: 12-10-2020
    """

    # load data
    data_array_int = np.fromfile(filename, np.dtype(">i"))
    data_array_float = np.fromfile(filename, np.dtype(">f"))

    # check version
    ver = data_array_int[0]
    if not ver == -1:
        sys.exit("incorrect version # "+str(ver)+" (not -1) found in file")

    # size of data array
    data_size = data_array_int[1]

    # reshape data into array
    data_array_int = data_array_int[2:]
    data_array_float = data_array_float[2:]
    data_array_int = data_array_int.reshape(data_size,4)
    data_array_float = data_array_float.reshape(data_size,4)

    # get vertex indices and coordinates
    ind = data_array_int[:,0]
    ind[ind < 0] = -ind[ind < 0] - 1
    ind[ind >= 0] = ind[ind >= 0] - 1

    x = data_array_float[:,1]
    y = data_array_float[:,2]
    z = data_array_float[:,3]

    return x, y, z, ind
