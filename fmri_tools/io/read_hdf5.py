# -*- coding: utf-8 -*-

# external inputs
import h5py
from nibabel.freesurfer.mghformat import MGHHeader


def read_hdf5(file_in):
    """Read HDF5.

    This function reads an hdf5 file which is expected to contain the datasets
    array, affine and header.

    Parameters
    ----------
    file_in : str
        File name of input file.

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
        Affine transformation matrix.
    header : MGHHeader
        Image header.

    """

    # check filename
    if isinstance(file_in, str):
        if not (file_in.endswith("h5") or file_in.endswith("hdf5")):
            raise ValueError("Currently supported file formats are " +
                             "h5 and hdf5.")
    else:
        raise ValueError("Filename must be a string!")

    with h5py.File(file_in, "r") as hf:

        # read data array
        if "array" in hf.keys():
            data = hf["array"][:]
        else:
            raise ValueError("No dataset found with name array!")

        # read affine matrix
        if "affine" in hf.keys():
            affine = hf["affine"][:]
        else:
            affine = None

        # read header
        if "header" in hf.keys():
            header = MGHHeader()
            header["dims"] = hf["header"]["dims"][:]
            header["Mdc"] = hf["header"]["Mdc"][:]
            header["Pxyz_c"] = hf["header"]["Pxyz_c"][:]
        else:
            header = None

    return data, affine, header
