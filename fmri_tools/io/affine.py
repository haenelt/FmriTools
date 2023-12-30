# -*- coding: utf-8 -*-
"""Header transformations."""

from pathlib import Path

import nibabel as nb
import numpy as np
from numpy.linalg import inv

__all__ = ["read_vox2vox", "read_vox2ras_tkr", "vox2ras_tkr"]


def read_vox2vox(input_lta):
    """This function reads a freesurfer lta file and extracts the vox2vox transformation
    matrix as numpy array.

    Parameters
    ----------
    input_lta : str
        Freesurfer lta file.

    Returns
    -------
    trans_formward : ndarray
        forwards affine transformation matrix.
    trans_inverse : ndarray
        inverse of M.

    """
    with open(input_lta, "r", encoding="utf-8") as f:
        x = f.readlines()

        # it is assumed that the vox2vox transformation matrix is found at
        # specific lines in the lta file
        transformation = x[8:12]

    # convert matrix to numpy array
    trans_forward = np.zeros((4, 4))
    for i, trans_ in enumerate(transformation):
        trans_forward[i, 0] = float(trans_.split()[0])
        trans_forward[i, 1] = float(trans_.split()[1])
        trans_forward[i, 2] = float(trans_.split()[2])
        trans_forward[i, 3] = float(trans_.split()[3])

        # get inverted transformation matrix (pseudoinverse)
        trans_inverse = np.linalg.pinv(trans_forward)

    return trans_forward, trans_inverse


def read_vox2ras_tkr(input_vol):
    """This function computes the vox2ras-tkr transform from header information of a
    Nifti or MGH file.

    Parameters
    ----------
    input_vol : str
        File name of Nifti or MGH volume.

    Returns
    -------
    ndarray
        Forward vox2ras-tkr transformation matrix.
    ndarray
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

    return vox2ras_tkr(dims, ds)


def vox2ras_tkr(dims, ds):
    """This function computes the vox2ras-tkr transform. The calculation follows the
    implementation found in the MGHHeader class
    (nibabel.freesurfer.mghformat.MGHHeader.get_vox2ras_tkr).

    Parameters
    ----------
    dims : tuple
        Tuple containing volume dimensions in x-, y- and z-direction.
    ds : tuple
        Tuple containing voxel sizes in x-, y- and z-direction.

    Returns
    -------
    v2rtkr : ndarray
        Forward vox2ras-tkr transformation matrix.
    r2vtkr : ndarray
        Inverse of v2rtkr.

    """
    ns = np.asarray(dims) * np.asarray(ds) / 2.0
    v2rtkr = np.array(
        [
            [-ds[0], 0, 0, ns[0]],
            [0, 0, ds[2], -ns[2]],
            [0, -ds[1], 0, ns[1]],
            [0, 0, 0, 1],
        ],
        dtype=np.float32,
    )
    r2vtkr = inv(v2rtkr)

    return v2rtkr, r2vtkr
