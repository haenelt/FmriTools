# -*- coding: utf-8 -*-
"""Apply transformations."""

import os
import subprocess

import nibabel as nb
import numpy as np
import numpy.linalg as npl
from numpy.matlib import repmat

from ..io.filename import get_filename

__all__ = ["apply_affine_chunked", "scanner_transform", "apply_warp", "apply_header"]


def apply_affine_chunked(aff, pts, chunk_size=10000):
    """Apply an affine matrix to points in chunks.

    This function is a copy of the routine `apply_affine` from the `affines` module of
    the `nibabel` package and applies an affine matrix to points in chunks. The only
    difference is that this function applies the affine transformation in chunks to
    prevent memory errors when working with large arrays. More information about this
    function can be found in the docstring of the aforementioned nibabel function.

    Parameters
    ----------
    aff : (N, N) np.ndarray
        Homogenous affine, for 3D points, will be 4 by 4. Contrary to first
        appearance, the affine will be applied on the left of `pts`.
    pts : (..., N-1) np.ndarray
        Points, where the last dimension contains the coordinates of each
        point. For 3D, the last dimension will be length 3.
    chunk_size : int, optional
        Chunk size for large arrays.

    Returns
    -------
    (..., N-1) np.ndarray
        Transformed points.

    """
    aff = np.asarray(aff)
    pts = np.asarray(pts)
    shape = pts.shape
    pts = pts.reshape((-1, shape[-1]))
    # rzs == rotations, zooms, shears
    rzs = aff[:-1, :-1]
    trans = aff[:-1, -1]

    # chunk intervals
    chunk = np.arange(0, len(pts), chunk_size)
    chunk[-1] = len(pts)
    if chunk[0] == 0:
        chunk = np.delete(chunk, 0)

    # apply affine transformation piecewise
    j1 = 0
    res = np.zeros_like(pts)
    for _, j2 in enumerate(chunk):
        res[j1:j2, :] = np.dot(pts[j1:j2, :], rzs.T) + trans[None, :]
        j1 = j2

    return res.reshape(shape)


def scanner_transform(input_source, input_target, path_output, compress_file=False):
    """This function uses the scanner coordinates to create a coordinate map between
    two images in the same scanner coordinate system. The orientation matrices written
    in the header of both files are taken to get the trasformation between both images.
    The output contains a 4d coordinate map describing the transformation from source to
    target image. Input files should be in nifti format.

    Parameters
    ----------
    input_source : str
        Absolute path of source file.
    input_target : str
        Absolute path of target file.
    path_output : str
        Path where output is saved.
    compress_file : bool, optional
        GZIP output. The default is False.

    Returns
    -------
    None.

    """
    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # load source and target images
    source_img = nb.load(input_source)
    target_img = nb.load(input_target)

    # get affine transformation
    target2source = npl.inv(source_img.affine).dot(target_img.affine)

    # initialise target coordinate map
    x_size = target_img.header["dim"][1]
    y_size = target_img.header["dim"][2]
    z_size = target_img.header["dim"][3]

    coordinate_mapping = np.zeros((x_size, y_size, z_size, 3), dtype="float")

    # coordinate mapping in x-direction
    X = np.array(np.arange(0, x_size, 1), dtype="float")
    X = np.transpose(repmat(X, y_size, 1))
    X = np.dstack([X] * z_size)

    # coordinate mapping in y-direction
    Y = np.array(np.arange(0, y_size), dtype="float")
    Y = repmat(Y, x_size, 1)
    Y = np.dstack([Y] * z_size)

    # coordinate mapping in z-direction
    Z = np.ones((x_size, y_size, z_size))
    Z = np.arange(0, z_size) * Z

    # merge directions
    coordinate_mapping[:, :, :, 0] = X
    coordinate_mapping[:, :, :, 1] = Y
    coordinate_mapping[:, :, :, 2] = Z

    # apply transformation to coordinates
    coordinate_mapping = apply_affine_chunked(target2source, coordinate_mapping)

    # write coordinate map
    target_img.header["dim"][0] = 4
    target_img.header["dim"][4] = 3
    target_img.set_data_dtype(np.float)

    # get filenames
    if os.path.splitext(os.path.basename(input_source))[1] == ".gz":
        name_source = os.path.splitext(
            os.path.splitext(os.path.basename(input_source))[0]
        )[0]
    else:
        name_source = os.path.splitext(os.path.basename(input_source))[0]

    if os.path.splitext(os.path.basename(input_target))[1] == ".gz":
        name_target = os.path.splitext(
            os.path.splitext(os.path.basename(input_target))[0]
        )[0]
    else:
        name_target = os.path.splitext(os.path.basename(input_target))[0]

    output = nb.Nifti1Image(coordinate_mapping, target_img.affine, target_img.header)
    if compress_file == True:
        nb.save(
            output,
            os.path.join(
                path_output, name_source + "_2_" + name_target + "_scanner.nii.gz"
            ),
        )
    else:
        nb.save(
            output,
            os.path.join(
                path_output, name_source + "_2_" + name_target + "_scanner.nii"
            ),
        )


def apply_warp(file_in, file_field, file_out):
    """Applies an FSL deformation field to a nifti image.

    Parameters
    ----------
    file_in : str
        File name of input image.
    file_field : str
        File name of deformation field.
    file_out : str
        File name of deformed image.
    """
    command = "applywarp"
    command += " --in=" + str(file_in)
    command += " --ref=" + str(file_in)
    command += " --out=" + str(file_out)
    command += " --warp=" + str(file_field)
    command += " --interp=spline --rel"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")


def apply_header(file_source, file_target, file_out, interp_method="nearest"):
    """Apply transformation to target image based on header information using
    FreeSurfer.

    Parameters
    ----------
    file_source : str
        File name of input image.
    file_target : str
        File name of target image.
    file_out : str
        File name of output image.
    interp_method : str, optional
        Interpolation method (nearest, trilin, cubic), by default "nearest"
    """
    # get filename
    path_out, _, _ = get_filename(file_out)

    # make output folder
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    command = "mri_vol2vol"
    command += " --no-save-reg"
    command += f" --interp {interp_method}"
    command += " --regheader"
    command += f" --mov {file_source}"
    command += f" --targ {file_target}"
    command += f" --o {file_out}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")
