# -*- coding: utf-8 -*-
"""Apply transformations."""

import os
import subprocess

import nibabel as nb
import numpy as np
import numpy.linalg as npl
from numpy.matlib import repmat

from ..io.filename import get_filename
from ..utils.interpolation import linear_interpolation3d, nn_interpolation3d

__all__ = [
    "apply_coordinate_mapping",
    "apply_affine_chunked",
    "scanner_transform",
    "apply_warp",
    "apply_flirt",
    "apply_header",
    "apply_fugue",
]


# linear and nearest neighbor sampling functions
_sampler = {"linear": linear_interpolation3d, "nearest": nn_interpolation3d}


def apply_coordinate_mapping(file_in, cmap_in, file_out, interpolation="linear"):
    """This function applies a coordinate mapping to a volume.

    Parameters
    ----------
    file_in : str
        Filename of input volume.
    cmap_in : str
        Filename of coordinate mapping.
    file_out : str
        Filename of output volume.
    interpolation : str, optional
        Interpolation type (linear or nearest). The default is "linear".

    Returns
    -------
    niimg
        Transformed volume.

    """
    # make output folder
    path_output = os.path.dirname(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # load data
    data = nb.load(file_in)
    arr = data.get_fdata()
    cmap = nb.load(cmap_in)
    arr_c = cmap.get_fdata()

    # get source and target image dimensions
    x_dim_source = data.header["dim"][1]
    y_dim_source = data.header["dim"][2]
    z_dim_source = data.header["dim"][3]
    t_dim_source = data.header["dim"][4]
    x_dim_target = cmap.header["dim"][1]
    y_dim_target = cmap.header["dim"][2]
    z_dim_target = cmap.header["dim"][3]

    # get mapping coordinates
    arr_c_x = arr_c[:, :, :, 0]
    arr_c_y = arr_c[:, :, :, 1]
    arr_c_z = arr_c[:, :, :, 2]

    # flatten
    arr_c_x = arr_c_x.flatten()
    arr_c_y = arr_c_y.flatten()
    arr_c_z = arr_c_z.flatten()

    # remove boundary coordinates
    arr_c_x = _set_min(arr_c_x, min_val=0.0)
    arr_c_y = _set_min(arr_c_y, min_val=0.0)
    arr_c_z = _set_min(arr_c_z, min_val=0.0)
    arr_c_x = _set_max(arr_c_x, max_val=x_dim_source - 1)
    arr_c_y = _set_max(arr_c_y, max_val=y_dim_source - 1)
    arr_c_z = _set_max(arr_c_z, max_val=z_dim_source - 1)

    # get coordinates to keep
    arr_sum = arr_c_x + arr_c_y + arr_c_z
    ind_ignore = np.where(np.isnan(arr_sum))[0]
    ind_keep = np.arange(len(arr_c_x))
    ind_keep = np.array(list(set(ind_keep).difference(set(ind_ignore))))

    # only use kept coordinates for interpolation
    arr_c_x = arr_c_x[ind_keep]
    arr_c_y = arr_c_y[ind_keep]
    arr_c_z = arr_c_z[ind_keep]

    # get dummy fourth dimension for 3d volumes to have support for 3D and 4D volumes.
    if t_dim_source == 1:
        arr = np.expand_dims(arr, axis=-1)

    # do the interpolation
    res = np.zeros((x_dim_target, y_dim_target, z_dim_target, t_dim_source))
    for i in range(t_dim_source):
        arr_sampled = _sampler[interpolation](
            arr_c_x, arr_c_y, arr_c_z, arr[:, :, :, i]
        )

        # reshape to output array
        tmp = np.zeros_like(arr_sum)
        tmp[ind_keep] = arr_sampled
        res[:, :, :, i] = np.reshape(tmp, (x_dim_target, y_dim_target, z_dim_target))

    # remove dummy dimensions in the case of 3D volume
    res = np.squeeze(res)

    output = nb.Nifti1Image(res, cmap.affine, cmap.header)
    nb.save(output, file_out)

    return output


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
    X = np.array(np.arange(0, x_size, 1), dtype=np.float64)
    X = np.transpose(repmat(X, y_size, 1))
    X = np.dstack([X] * z_size)

    # coordinate mapping in y-direction
    Y = np.array(np.arange(0, y_size), dtype=np.float64)
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
    target_img.set_data_dtype(np.float64)

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


def apply_flirt(file_in, file_ref, file_mat, file_out, interp_method="trilinear"):
    """Apply transformation from flirt registration to image.

    Parameters
    ----------
    file_in : str
        File name of input file.
    file_ref : str
        File name of reference file.
    file_mat : str
        File name of 4x4 affine matrix saved as *.txt or *.mat file.
    file_out : str
        File name of transformed input file.
    interp_method : str, optional
        Interpolation method (trilinear, nearestneighbor, sinc, spline), by default
        "trilinear"
    """
    # get filename
    path_out, _, _ = get_filename(file_out)

    # make output folder
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    command = "flirt"
    command += f" -in {file_in}"
    command += f" -ref {file_ref}"
    command += f" -out {file_out}"
    command += " -omat inn_flirt.mat -applyxfm"
    command += f" -init {file_mat}"
    command += f" -interp {interp_method}"
    command += " -paddingsize 0"

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


def apply_fugue(file_in, file_shift, udir, forward_warping=False):
    """Apply field map deformation to image using FSL.

    Parameters
    ----------
    file_in : _type_
        File name of input file.
    file_shift : _type_
        File name of shift file containing deformation.
    udir : _type_
        Unwarping direction (x, y, z, x-, y-, or z-).
    forward_warping : bool, optional
        Apply forward warping instead of unwarping, by default False
    """
    # output file name
    _, name_in, ext_in = get_filename(file_in)

    command = "fugue"
    command += f" --in={file_in}"
    command += f" --loadshift={file_shift}"
    command += f" --unwarpdir={udir}"
    if forward_warping:
        command += f" --warp={name_in}_warped{ext_in}"
    else:
        command += f" --unwarp={name_in}_unwarped{ext_in}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")


def _set_min(arr, min_val):
    """Remove coordinates below the matrix size."""
    arr[np.floor(arr) < min_val] = None
    return arr


def _set_max(arr, max_val):
    """Remove coordinates above the matrix size."""
    arr[np.ceil(arr) > max_val] = None
    return arr
