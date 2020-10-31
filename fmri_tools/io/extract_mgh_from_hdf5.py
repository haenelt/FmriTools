# -*- coding: utf-8 -*-

# external inputs
import numpy as np
from nibabel.freesurfer.mghformat import MGHHeader

# local input
from fmri_tools.io import read_hdf5
from fmri_tools.io import write_mgh


def extract_mgh_from_hdf5(file_in, file_out, t, n):
    """ Extract MGH from HDF5
    
    This function reads an hdf5 file which is expected to contain a 3D array 
    with dimensions vertex x time point x layer and optionally some header
    information and an affine transformation matrix from the original mgh file. 
    Data from one time point and one layer are extracted and saved again as mgh 
    file. If no affine matrix or header information exist, an identity matrix 
    and an empty header are set, respectively.

    Parameters
    ----------
    file_in : str
        Filename of hdf5 input file.
    file_out : str
        Filename of mgh output file.
    t : int
        Time point which is assumed to be stored along the second dimension.
    n : int
        Layer which is assumed to be stored along the third dimension.

    Raises
    ------
    ValueError
        If the read 'data' array has not the right number of dimensions.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 20-10-2020
    Last modified: 31-10-2020

    """    
    
    # read file
    data, affine, header = read_hdf5(file_in)
    
    # check dimensionality
    if len(np.shape(data)) != 3:
        raise ValueError("Data array has incorrect number of dimensions!")
        
    # extract one time point and one layer
    data = data[:,t,n]
        
    # check affine
    if affine is None:
        affine = np.eye(4)
    
    # check header
    if header is None:
        header = MGHHeader()
    
    # write MGH file
    write_mgh(file_out, data, affine, header)
