# -*- coding: utf-8 -*-
"""Input/Output surface utilities."""

import datetime
import os
import shutil as sh
import subprocess
import sys
from pathlib import Path

import nibabel as nb
import numpy as np
from nibabel.freesurfer.io import (
    read_label,
    read_morph_data,
    write_geometry,
    write_morph_data,
)
from nibabel.freesurfer.mghformat import MGHHeader
from scipy.spatial import Delaunay

from .filename import get_filename

__all__ = [
    "write_mgh",
    "read_mgh",
    "write_label",
    "read_patch",
    "patch_as_mesh",
    "mgh_to_patch",
    "curv_to_patch",
    "label_to_patch",
    "label_as_patch",
]


def write_mgh(file_out, arr, affine=None, header=None):
    """This function adds two empty dimensions to an array and saves it as a freesurfer
    mgh surface file.

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

    if not str(file_out).endswith("mgh"):
        raise ValueError("Currently supported file format is mgh.")

    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # add empty dimensions
    arr = np.expand_dims(arr, axis=1)
    arr = np.expand_dims(arr, axis=1)

    if affine is None:
        affine = np.eye(4)

    if header is None:
        header = MGHHeader()

    # write output
    output = nb.Nifti1Image(arr, affine, header)
    nb.save(output, file_out)


def read_mgh(file_in):
    """This function reads a surface mgh file and removes empty dimensions from the data
    array.

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

    if not str(file_in).endswith("mgh"):
        raise ValueError("Currently supported file format is mgh.")

    # get header
    header = nb.load(file_in).header
    affine = nb.load(file_in).affine

    # get data
    arr = nb.load(file_in).get_fdata()
    arr = np.squeeze(arr)

    return arr, affine, header


def write_label(file_out, arr_label):
    """This function writes a textfile which can be read as label file in freesurfer.

    Parameters
    ----------
    file_out : str
        Filename of label file.
    arr_label : list
        List of label indices.

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

    if not str(file_out).endswith("label"):
        raise ValueError("Currently supported file format is txt.")

    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # number of labels
    n_label = len(arr_label)

    with open(file_out, "w", encoding="utf-8") as f:
        f.write("#!ascii label  , from subject  vox2ras=TkReg\n")
        f.write(str(n_label) + "\n")
        for i in range(n_label):
            f.write(str(arr_label[i]) + " 0.000 0.000 0.000 0.000\n")


def read_patch(file_in):
    """This function reads an freesurfer patch saved in binary format. It is an
    equivalent to the matlab function read_patch in the ./freesurfer/matlab folder. Data
    is read in big endian order.

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
    if ver != -1:
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


def patch_as_mesh(file_out, file_patch):
    """This function reads coordinates from a flattened freesurfer patch and creates a
    triangular surface mesh using Delaunay triangulation.

    Parameters
    ----------
    file_out : str
        Filename of output surface mesh.
    file_patch : str
        Filename of freesurfer patch file.

    Returns
    -------
    None.

    """
    x, y, _, _ = read_patch(file_patch)
    coords = np.zeros((len(x), 2))
    coords[:, 0] = x
    coords[:, 1] = y

    tri = Delaunay(coords)
    fac = tri.simplices

    vtx = np.zeros((len(x), 3))
    vtx[:, 0] = x
    vtx[:, 1] = y

    write_geometry(file_out, vtx, fac)


def mgh_to_patch(file_out, file_mgh, file_patch):
    """This function reads an MGH overlay and saves only coordinates that match with the
    corresponding freesurfer patch.

    Parameters
    ----------
    file_out : str
        Filename of output surface mesh.
    file_mgh : str
        Filename of mgh overlay.
    file_patch : str
        Filename of freesurfer patch file.

    Returns
    -------
    None.

    """
    _, _, _, ind = read_patch(file_patch)
    arr, affine, header = read_mgh(file_mgh)
    arr = arr[ind]
    write_mgh(file_out, arr, affine=affine, header=header)


def curv_to_patch(file_out, file_curv, file_patch):
    """This function reads a freesurfer curvature file and saves only coordinates that
    match with the corresponding freesurfer patch.

    Parameters
    ----------
    file_out : str
        Filename of output surface mesh.
    file_curv : str
        Filename of freesurfer curvature file.
    file_patch : str
        Filename of freesurfer patch file.

    Returns
    -------
    None.

    """
    _, _, _, ind = read_patch(file_patch)
    curv = read_morph_data(file_curv)
    curv = curv[ind]
    write_morph_data(file_out, curv)


def label_to_patch(file_out, file_label, file_patch):
    """This function reads a freesurfer label file and saves only coordinates that match
    with the corresponding freesurfer patch.

    Parameters
    ----------
    file_out : str
        Filename of output surface mesh.
    file_label : str
        Filename of freesurfer label file.
    file_patch : str
        Filename of freesurfer patch file.

    Returns
    -------
    None.

    """
    _, _, _, ind = read_patch(file_patch)
    label = read_label(file_label)
    label_new = [np.where(ind == i)[0][0] for i in label]
    write_label(file_out, label_new)


def label_as_patch(file_ref, file_label, file_out, cleanup=True):
    """Uses the FreeSurfer label2patch function to convert a label file into a *.path.3d
    file. This can be useful for surface flattening of a region defined by a label file.

    Parameters
    ----------
    file_ref : str
        Reference surface file.
    file_label : str
        Label file.
    file_output : str
        File name of patch file.
    cleanup : bool, optional
        Delete intermediate files. The default is True.

    Returns
    -------
    None.

    """
    # create temporary folder
    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = "".join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    path_temp = os.path.join(os.path.dirname(file_out), "tmp_" + tmp_string)
    path_subj = os.path.join(path_temp, "subj")
    path_surf = os.path.join(path_subj, "surf")

    # make temporary folder
    if not os.path.exists(path_temp):
        os.makedirs(path_temp)
        os.mkdir(path_subj)
        os.mkdir(path_surf)
    else:
        raise FileExistsError("Temporary folder already exists!")

    # get hemisphere from file name
    hemi = os.path.basename(file_label)[:2]

    # copy reference file into temporary folder
    sh.copy2(file_ref, os.path.join(path_surf, hemi + "." + "surf"))

    # convert label to patch
    command = "label2patch"
    command += " -sdir " + str(path_temp)
    command += " -surf surf"
    command += " subj"
    command += " " + str(hemi)
    command += " " + str(file_label)
    command += " " + str(file_out)

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")

    # delete temporary files
    if cleanup:
        sh.rmtree(path_temp, ignore_errors=True)
