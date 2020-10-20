# -*- coding: utf-8 -*-

# external inputs
import numpy as np
from nibabel.freesurfer.mghformat import MGHHeader

# local input
from fmri_tools.io import read_hdf5
from fmri_tools.io import write_mgh


def extract_mgh_from_hdf5(file_in, file_out, t=0, n=0):
    """ Extract MGH from HDF5
    
    This function reads an hdf5 file which is expected to contain a 3D array 
    containing data points vertex x time point x layer and extract data for
    one time point and one layer. The extracted data is saved as mgh file.

    Parameters
    ----------
    file_in : str
        Filename of hdf5 input file.
    file_out : str
        Filename of mgh output file.
    t : int, optional
        Time point which is assumed to be stored along the second dimension. The 
        default is 0.
    n : int, optional
        Layer which is assumed to be stored along the third dimension. The 
        default is 0.

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
    Last modified: 20-10-2020

    """    
    
    # read file
    data, affine, header = read_hdf5(file_in)
    
    # check dimensionality
    if len(np.shape(data)) != 3:
        raise ValueError("Data array has incorrect number of dimensions!")
        
    # extract one time point and one layer
    data = data[:,t,n]
        
    # check affine
    if not affine:
        affine = np.eye(4)
    
    # check header
    if not header:
        header = MGHHeader()
    
    # write MGH file
    write_mgh(file_out, data, affine, header)
