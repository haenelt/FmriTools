# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import h5py
import numpy as np
from nibabel.freesurfer.mghformat import MGHHeader

# local input
from fmri_tools.io.get_filename import get_filename


def write_hdf5(file_out, arr, affine=None, header=None):
    """ Write HDF5
    
    This function writes a numpy array to an hdf5 file. Optionally, an affine 
    transformation matrix and parts of an MGHHeader are stored. The array is 
    stored in half-precision floating-point format.

    Parameters
    ----------
    file_out : str
        Filename of output file.
    arr : ndarray
        Image array.
    affine : ndarray, optional
        Affine transformation matrix. The default is None.
    header : MGHHeader, optional
        Image header. The default is None.

    Raises
    ------
    ValueError
        If `file_out` is not a string or has a file extension which is not 
        supported.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 20-10-2020
    Last modified: 25-10-2020

    """

    # check filename
    if isinstance(file_out, str):
        if not (file_out.endswith("h5") or file_out.endswith("hdf5")):            
            raise ValueError("Currently supported file formats are " + \
                             "h5 and hdf5.")
    else:
        raise ValueError("Filename must be a string!")

    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # create new MGH header
    if header is not None:
        header_new = MGHHeader()
        header_new["dims"][0] = len(arr)
        header_new["Mdc"] = header["Mdc"]
        header_new["Pxyz_c"] = header["Pxyz_c"]

    with h5py.File(file_out, "w") as hf:
        hf.create_dataset("array",  
                          data=arr, 
                          shape=np.shape(arr),
                          chunks=True,
                          compression="gzip", 
                          compression_opts=9,
                          dtype=np.float16)
        
        if affine is not None:
            hf.create_dataset("affine", 
                              data=affine,
                              shape=np.shape(affine),
                              compression="gzip", 
                              compression_opts=9)

        if header is not None:
            grp = hf.create_group("header")
            grp.create_dataset("dims", 
                               data=header_new["dims"], 
                               compression="gzip", 
                               compression_opts=9)
            grp.create_dataset("Mdc", 
                               data=header_new["Mdc"], 
                               compression="gzip", 
                               compression_opts=9)
            grp.create_dataset("Pxyz_c", 
                               data=header_new["Pxyz_c"], 
                               compression="gzip", 
                               compression_opts=9)
