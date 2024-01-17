# -*- coding: utf-8 -*-
"""FSL registration algorithms."""

import os
import shutil as sh
import subprocess

from ..io.filename import get_filename
from ..io.vol import mri_convert
from ..registration.cmap import generate_coordinate_mapping
from ..registration.transform import (
    apply_coordinate_mapping,
    apply_flirt,
    scanner_transform,
)

__all__ = ["flirt", "get_flash2orig"]


def flirt(
    file_in,
    file_ref,
    file_out,
    file_out_matrix,
    cost_func="mutualinfo",
    interp_method="trilinear",
):
    """Run FSL affine flirt registration.

    Parameters
    ----------
    file_in : str
        File name of input image (source).
    file_ref : str
        File name of reference image (target).
    file_out : str
        File name of transformed output image.
    file_out_matrix : str
        File name of text file (*.txt) containing the estimated transformation matrix.
    cost_func : str, optional
        Cost function for regiration ('mutualinfo', 'corratio', 'normcorr', 'normmi',
        'leastsq', 'labeldiff' or 'bbr')), by default "mutualinfo"
    interp_method : str, optional
        Interpolation method, by default "trilinear"

    Raises
    ------
    ValueError
        If file_out_matrix does not have an *.txt extension.
    """
    # make output folder
    path_out, _, _ = get_filename(file_out)
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    # check matrix is txt
    _, _, ext_matrix = get_filename(file_out_matrix)
    if ext_matrix != ".txt":
        raise ValueError("file_out_matrix must have a text file (*.txt) extension!")

    command = "flirt"
    command += f" -in {file_in}"
    command += f" -ref {file_ref}"
    command += f" -out {file_out}"
    command += f" -omat {file_out_matrix}"
    command += f" -searchcost {cost_func}"
    command += " -dof 6"
    command += f" -interp {interp_method}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")


def get_flash2orig(file_flash, file_inv2, file_orig, path_output, cleanup=False):
    """This function computes the deformation field for the registration between a
    partial coverage GRE image and the freesurfer orig file. The following steps are
    performed: (1) set output folder structure, (2) get scanner transform GRE <-> inv2
    and inv2 -> orig, (3) generate flash cmap, (4) apply scanner transform inv2 -> GRE,
    (5) get flirt registration GRE -> inv2, (6) apply flirt to GRE cmap, (7) apply
    scanner transform to GRE cmap, (8) apply final deformation to GRE. The function
    needs the FSL environment set.

    Parameters
    ----------
    file_flash : str
        Input path for GRE image.
    file_inv2 : str
        Input path for MP2RAGE INV2 image.
    file_orig : str
        Input path for freesurfer orig image.
    path_output : str
        Path where output is saved.
    cleanup : bool, optional
        Delete intermediate files. The default is False.

    Returns
    -------
    None.

    """
    # set folder structure
    path_temp = os.path.join(path_output, "temp")

    if not os.path.exists(path_output):
        os.makedirs(path_output)

    if not os.path.exists(path_temp):
        os.makedirs(path_temp)

    # copy input files
    sh.copyfile(file_inv2, os.path.join(path_temp, "inv2.nii"))
    sh.copyfile(file_flash, os.path.join(path_temp, "flash.nii"))

    # convert orig to nifti
    mri_convert(file_orig, os.path.join(path_temp, "orig.nii"))

    # scanner transformation
    scanner_transform(
        os.path.join(path_temp, "inv2.nii"),
        os.path.join(path_temp, "flash.nii"),
        path_temp,
        False,
    )
    scanner_transform(
        os.path.join(path_temp, "flash.nii"),
        os.path.join(path_temp, "inv2.nii"),
        path_temp,
        False,
    )
    scanner_transform(
        os.path.join(path_temp, "inv2.nii"),
        os.path.join(path_temp, "orig.nii"),
        path_temp,
        False,
    )

    # generate coordinate mapping
    generate_coordinate_mapping(
        os.path.join(path_temp, "flash.nii"), 0, path_temp, "flash", False, True
    )

    # scanner transform inv2 to flash
    apply_coordinate_mapping(
        os.path.join(path_temp, "inv2.nii"),
        os.path.join(path_temp, "inv2_2_flash_scanner.nii"),
        os.path.join(path_temp, "inv2_apply_scanner_def-img.nii.gz"),
        interpolation="linear",
    )

    # flirt flash to inv2
    os.chdir(path_temp)
    flirt(
        os.path.join(path_temp, "flash.nii"),
        os.path.join(path_temp, "inv2_apply_scanner_def-img.nii.gz"),
        os.path.join(path_temp, "flash_apply_flirt_def-img.nii.gz"),
        os.path.join(path_temp, "flirt_matrix.txt"),
        cost_func="mutualinfo",
        interp_method="trilinear",
    )

    # apply flirt to flash cmap
    apply_flirt(
        os.path.join(path_temp, "cmap_flash.nii"),
        os.path.join(path_temp, "flash.nii"),
        os.path.join(path_temp, "flirt_matrix.txt"),
        os.path.join(path_temp, "cmap_apply_flirt_def-img.nii.gz"),
        "trilinear",
    )

    # concatenate cmaps
    apply_coordinate_mapping(
        os.path.join(path_temp, "flash_2_inv2_scanner.nii"),
        os.path.join(path_temp, "inv2_2_orig_scanner.nii"),
        os.path.join(path_temp, "flash_2_orig_scanner.nii"),
        interpolation="linear",
    )

    # apply scanner transform to flash cmap
    apply_coordinate_mapping(
        os.path.join(path_temp, "cmap_apply_flirt_def-img.nii.gz"),
        os.path.join(path_temp, "flash_2_orig_scanner.nii"),
        os.path.join(path_temp, "cmap_apply_scanner_def-img.nii.gz"),
        interpolation="linear",
    )

    # apply deformation to source image
    apply_coordinate_mapping(
        os.path.join(path_temp, "flash.nii"),  # input
        os.path.join(path_temp, "cmap_apply_scanner_def-img.nii.gz"),
        os.path.join(path_temp, "flash_apply_deformation_def-img.nii.gz"),
        interpolation="linear",  # nearest or linear
    )

    # rename final deformation examples
    os.rename(
        os.path.join(path_temp, "cmap_apply_scanner_def-img.nii.gz"),
        os.path.join(path_output, "flash2orig.nii.gz"),
    )
    os.rename(
        os.path.join(path_temp, "flash_apply_deformation_def-img.nii.gz"),
        os.path.join(path_output, "flash2orig_example.nii.gz"),
    )

    # clean intermediate files
    if cleanup:
        sh.rmtree(path_temp, ignore_errors=True)
