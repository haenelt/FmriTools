# -*- coding: utf-8 -*-

# python standard library inputs
import os
from pathlib import Path

# external inputs
import numpy as np
import h5py
from nibabel.freesurfer.mghformat import MGHHeader

# local input
from .surf import write_mgh
from .get_filename import get_filename

__all__ = ['read_hdf5', 'write_hdf5', 'extract_mgh_from_hdf5']


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
    if not isinstance(file_in, str) and not isinstance(file_in, Path):
        raise ValueError("Filename must be a string or a pathlib.Path instance!")

    if not str(file_in).endswith("h5") and not str(file_in).endswith("hdf5"):
        raise ValueError("Currently supported file formats are " +
                         "h5 and hdf5.")

    with h5py.File(file_in, "r") as hf:
        # read data array
        if "array" in hf.keys():
            data = hf["array"][:]
        else:
            raise ValueError("No dataset found with name array!")

        # read affine matrix
        affine = hf["affine"][:] if "affine" in hf.keys() else None
        # read header
        if "header" in hf.keys():
            header = MGHHeader()
            header["dims"] = hf["header"]["dims"][:]
            header["Mdc"] = hf["header"]["Mdc"][:]
            header["Pxyz_c"] = hf["header"]["Pxyz_c"][:]
        else:
            header = None

    return data, affine, header


def write_hdf5(file_out, arr, affine=None, header=None):
    """Write HDF5.

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

    """

    # check filename
    if not isinstance(file_out, str) and not isinstance(file_out, Path):
        raise ValueError("Filename must be a string or a pathlib.Path instance!")

    if not str(file_out).endswith("h5") and not str(file_out).endswith("hdf5"):
        raise ValueError("Currently supported file formats are " +
                         "h5 and hdf5.")

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


def extract_mgh_from_hdf5(file_in, file_out, t, n):
    """Extract MGH from HDF5.

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

    """

    # read file
    data, affine, header = read_hdf5(file_in)

    # check dimensionality
    if len(np.shape(data)) != 3:
        raise ValueError("Data array has incorrect number of dimensions!")

    # extract one time point and one layer
    data = data[:, t, n]

    # check affine
    if affine is None:
        affine = np.eye(4)

    # check header
    if header is None:
        header = MGHHeader()

    # write MGH file
    write_mgh(file_out, data, affine, header)
