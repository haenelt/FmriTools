# -*- coding: utf-8 -*-
"""Gradient nonlinearity correction."""

import os
import shutil as sh
import subprocess

import nibabel as nb
import numpy as np

from ..io.filename import get_filename
from ..registration.cmap import generate_coordinate_mapping
from ..registration.transform import apply_warp, combine_warp, convert_warp
from ..utils.calc import mean_image


def gnl_correction(
    file_in, file_bash, file_coeff, python3_env, python2_env, path_output, cleanup=True
):
    """The purpose of the following function is to correct for gradient nonlinearities.
    A corrected file is written using spline interpolation. The function needs FSL to be
    included in the search path.

    Parameters
    ----------
    file_in : str
        Filename of input image.
    file_bash : str
        Filename of bash script which calls the gradient unwarping toolbox.
    file_coeff : str
        Filename of siemens coefficient file.
    python3_env : str
        Name of python3 virtual environment.
    python2_env : str
        Name of python2 virtual environment.
    path_output : str
        Path where output is written.
    cleanup : bool, optional
        Delete intermediate files. The default is True.

    Returns
    -------
    None.

    """
    # get fileparts
    _, name, ext = get_filename(file_in)

    # make subfolders
    path_grad = os.path.join(path_output, "grad")
    if not os.path.exists(path_grad):
        os.makedirs(path_grad)

    # parse arguments
    file_output = os.path.join(path_output, name + "_gnlcorr" + ext)
    file_warp = os.path.join(path_grad, "warp.nii.gz")
    file_jacobian = os.path.join(path_grad, "warp_jacobian.nii.gz")

    # run gradient unwarp
    command = "bash"
    command += f" {file_bash}"
    command += f" {python3_env}"
    command += f" {python2_env}"
    command += f" {path_grad}"
    command += f" {file_in}"
    command += " trilinear.nii.gz"
    command += f" {file_coeff}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")

    # now create an appropriate warpfield output (relative convention)
    convert_warp(
        os.path.join(path_grad, "trilinear.nii.gz"),
        os.path.join(path_grad, "fullWarp_abs.nii.gz"),
        file_jacobian,
        os.path.join(path_grad, "fullWarp_abs.nii.gz"),
    )

    # convertwarp's jacobian output has 8 frames, each combination of one-sided
    # differences, so average them
    mean_image(file_jacobian, file_jacobian)

    # apply warp to first volume
    apply_warp(file_in, file_warp, file_output)

    # normalise warped output image to initial intensity range
    data_img = nb.load(file_in)
    data_array = data_img.get_fdata()
    max_data = np.max(data_array)
    min_data = np.min(data_array)

    data_img = nb.load(file_output)
    data_array = data_img.get_fdata()
    data_array[data_array < min_data] = 0
    data_array[data_array > max_data] = max_data

    output = nb.Nifti1Image(data_array, data_img.affine, data_img.header)
    nb.save(output, file_output)

    # calculate gradient deviations
    command = "calc_grad_perc_dev"
    command += f" --fulwarp={file_warp}"
    command += f" -o {os.path.join(path_grad, 'grad_dev')}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")

    # merge directions
    combine_warp(
        os.path.join(path_grad, "grad_dev_x.nii.gz"),
        os.path.join(path_grad, "grad_dev_y.nii.gz"),
        os.path.join(path_grad, "grad_dev_z.nii.gz"),
        os.path.join(path_grad, "grad_dev.nii.gz"),
    )

    # convert from % deviation to absolute
    data_img = nb.load(os.path.join(path_grad, "grad_dev.nii.gz"))
    data_array = data_img.get_fdata()
    data_array = data_array / 100

    output = nb.Nifti1Image(data_array, data_img.affine, data_img.header)
    nb.save(output, os.path.join(path_grad, "grad_dev.nii.gz"))

    # warp coordinate mapping
    generate_coordinate_mapping(
        file_in, 0, path_grad, suffix="gnl", time=False, write_output=True
    )

    apply_warp(
        os.path.join(path_grad, "cmap_gnl.nii"),
        file_warp,
        os.path.join(path_grad, "cmap_gnl.nii"),
    )

    # clean intermediate files
    if cleanup:
        sh.rmtree(path_grad, ignore_errors=True)
