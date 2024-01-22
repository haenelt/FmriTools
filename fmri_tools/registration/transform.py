# -*- coding: utf-8 -*-
"""Apply transformations."""

import datetime
import os
import shutil as sh
import subprocess
from shutil import copyfile

import nibabel as nb
import numpy as np
import numpy.linalg as npl
from numpy.matlib import repmat
from sh import gunzip

from ..io.affine import apply_affine_chunked
from ..io.filename import get_filename
from ..utils.interpolation import linear_interpolation3d, nn_interpolation3d

__all__ = [
    "apply_coordinate_mapping",
    "scanner_transform",
    "apply_header",
    "resample_volume",
    "transform_timeseries",
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
        Interpolation method (nearest, trilin, cubic), by default "nearest".
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
        subprocess.run([command], shell=True, check=False)
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


def resample_volume(file_in, file_out, dxyz=(0.4, 0.4, 0.4), rmode="Cu"):
    """This function resamples a nifti volume using the afni function 3dresample. Before
    running the function, set the afni environment by calling AFNI in the terminal.

    Parameters
    ----------
    file_in : str
        Nifti input filename.
    file_out : str
        Nifti output filename.
    dxyz : tuple, optional
        Array of target resolution in single dimensions. The default is
        (0.4, 0.4, 0.4).
    rmode : str, optional
        Interpolation methods (Linear, NN, Cu, Bk). The default is "Cu".
    """
    # get path and file extension of input file
    path_in, _, ext_in = get_filename(file_in)

    # make temporary copy of input file
    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = "".join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    file_tmp = os.path.join(path_in, "tmp_" + tmp_string + ext_in)

    if not os.path.exists(file_tmp) and not os.path.exists(file_tmp[:-3]):
        copyfile(file_in, file_tmp)
    else:
        raise FileExistsError("Temporary file already exists!")

    if os.path.splitext(file_tmp)[1] == ".gz":
        gunzip(file_tmp)
        file_tmp = os.path.splitext(file_tmp)[0]

    # resample volume
    command = "3dresample"
    command += f" -dxyz {dxyz[0]} {dxyz[1]} {dxyz[2]}"
    command += f" -rmode {rmode}"
    command += f" -inset {file_tmp}"
    command += f" -prefix {file_out}"

    print("Execute: " + command)
    try:
        subprocess.run([command], shell=True, check=False)
    except subprocess.CalledProcessError:
        print("Execuation failed!")

    # remove temporary copy
    os.remove(file_tmp)


def transform_timeseries(file_in, file_cmap, interpolation):
    """Epi time series in native space are transformed to a target space using a
    deformation field. The transformed time series get the prefix r. Because of heap
    size limits, the time series is split and the deformation is applied separately to
    each volume.

    Parameters
    ----------
    file_in : str
        File name of time series.
    file_cmap : str
        File name of coordinate mapping.
    interpolation : str
        Interpolation method (linear or nearest).

    Raises
    ------
    FileExistsError
        If temporary folder already exists.
    """
    # make temporary output folder
    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = "".join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    path_tmp = os.path.join(os.path.dirname(file_in), "tmp_" + tmp_string)

    if not os.path.exists(path_tmp):
        os.mkdir(path_tmp)
    else:
        raise FileExistsError("Temporary folder already exists!")

    # get length of time series
    data = nb.load(file_in)
    data_array = data.get_fdata()
    nt = nb.load(file_in).header["dim"][4]
    data.header["dim"][4] = 1  # change header for single 3d volumes

    for j in range(nt):
        # save single time frame
        output = nb.Nifti1Image(data_array[:, :, :, j], data.affine, data.header)
        nb.save(output, os.path.join(path_tmp, f"{j}.nii"))

        apply_coordinate_mapping(
            os.path.join(path_tmp, f"{j}.nii"),
            file_cmap,
            os.path.join(path_tmp, f"{j}_def.nii"),
            interpolation=interpolation,
        )

    # merge final deformed time series
    data = nb.load(os.path.join(path_tmp, "0_def.nii"))
    data.header["dim"][4] = nt
    data_res = np.zeros(data.header["dim"][1:5])

    for j in range(nt):
        data_res[:, :, :, j] = nb.load(
            os.path.join(path_tmp, str(j) + "_def.nii")
        ).get_fdata()

        # time series path and basename
    path = os.path.dirname(file_in)
    file_ = os.path.splitext(os.path.basename(file_in))[0]

    output = nb.Nifti1Image(data_res, data.affine, data.header)
    nb.save(output, os.path.join(path, f"r{file_}_linear.nii"))

    # delete intermediate files
    sh.rmtree(path_tmp, ignore_errors=True)
