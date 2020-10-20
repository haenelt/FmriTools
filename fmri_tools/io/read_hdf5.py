# -*- coding: utf-8 -*-

# external inputs
import h5py
from nibabel.freesurfer.mghformat import MGHHeader


def read_hdf5(file_in):
    """ Read HDF5
    
    This function    

    Parameters
    ----------
    file_in : str
        Filename of input file.

    Raises
    ------
    ValueError
        If `file_in` is not a string or has a file extension which is not 
        supported.

    Returns
    -------
    arr : ndarray
        Image array.
    affine : ndarray
        Affine transformation matrix. The default is None.
    header : MGHHeader
        Image header. The default is None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 20-10-2020
    Last modified: 20-10-2020

    """
    
    # check filename
    if isinstance(file_in, str):
        if not (file_in.endswith("hdf") or file_in.endswith("h4") or
                file_in.endswith("hdf4") or file_in.endswith("he2") or 
                file_in.endswith("h5") or file_in.endswith("hdf5") or 
                file_in.endswith("he5")):            
            raise ValueError("Currently supported file formats are " + \
                             "hdf, h4, hdf4, he2, h5, hdf5 and he5.")
    else:
        raise ValueError("Filename must be a string!")

    with h5py.File(file_in, "r") as hf:
        
        # read data array
        if "array" in hf.keys():
            data = hf["array"][:]
        else:
            raise ValueError("No dataset with name array found!")
        
        # read affine matrix
        if "affine" in hf.keys():
            affine = hf["affine"][:]
        else:
            affine = []

        # read header
        if "header" in hf.keys():
            header = MGHHeader()
            header["dims"] = hf["header"]["dims"]
            header["Mdc"] = hf["header"]["Mdc"]
            header["Pxyz_c"] = hf["header"]["Pxyz_c"]
        else:
            header = []

    return data, affine, header
