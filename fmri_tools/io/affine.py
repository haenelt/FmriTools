# -*- coding: utf-8 -*-

# python standard library inputs
from pathlib import Path

# external inputs
import numpy as np
import nibabel as nb
from numpy.linalg import inv

__all__ = ['read_vox2vox', 'read_vox2ras_tkr']


def read_vox2vox(input_lta):
    """Read vox2vox.

    This function reads a freesurfer lta file and extracts the vox2vox
    transformation matrix as numpy array.

    Parameters
    ----------
    input_lta : str
        Freesurfer lta file.

    Returns
    -------
    M : ndarray
        forwards affine transformation matrix.
    Minv : ndarray
        inverse of M.

    """

    with open(input_lta, "r") as f:
        x = f.readlines()

        # it is assumed that the vox2vox transformation matrix is found at
        # specific lines in the lta file
        transformation = x[8:12]

    # convert matrix to numpy array
    M = np.zeros((4, 4))
    for i in range(len(transformation)):
        x = transformation[i].split()

        M[i, 0] = float(x[0])
        M[i, 1] = float(x[1])
        M[i, 2] = float(x[2])
        M[i, 3] = float(x[3])

        # get inverted transformation matrix (pseudoinverse)
        Minv = np.linalg.pinv(M)

    return M, Minv


def read_vox2ras_tkr(input_vol):
    """Calculate vox2ras_tkr from header information.

    This function computes the vox2ras-tkr transform from header information of a Nifti
    of MGH file. The calculation follows the implementation found in the MGHHeader class
    (nibabel.freesurfer.mghformat.MGHHeader.get_vox2ras_tkr).

    Parameters
    ----------
    input_vol : str
        File name of Nifti or MGH volume.

    Returns
    -------
    v2rtkr : ndarray
        Forward vox2ras-tkr transformation matrix.
    r2vtkr : ndarray
        Inverse of v2rtkr.

    """
    hdr = nb.load(input_vol).header
    ext = Path(input_vol).suffixes
    if ".nii" in ext:
        dims = hdr["dim"][1:4]
        ds = hdr["pixdim"][1:4]
    elif ".mgh" in ext or ".mgz" in ext:
        dims = hdr["dims"][:3]
        ds = hdr["delta"][:3]
    else:
        raise ValueError("Unknown file extension!")

    ns = dims * ds / 2.0
    v2rtkr = np.array(
        [[-ds[0], 0, 0, ns[0]],
         [0, 0, ds[2], -ns[2]],
         [0, -ds[1], 0, ns[1]],
         [0, 0, 0, 1]],
        dtype=np.float32,
    )
    r2vtkr = inv(v2rtkr)

    return v2rtkr, r2vtkr
