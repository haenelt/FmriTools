# -*- coding: utf-8 -*-
"""FSL registration algorithms."""

import os

from .. import execute_command
from ..io.filename import get_filename

__all__ = [
    "flirt",
    "apply_warp",
    "apply_flirt",
    "convert_warp",
    "combine_warp",
    "apply_fugue",
]


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

    # run
    execute_command(command)


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

    # run
    execute_command(command)


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

    # run
    execute_command(command)


def convert_warp(file_ref, file_warp, file_jacobian, file_out):
    """Conversion of warp field using using FSL to get Jacobian output.

    Parameters
    ----------
    file_ref : str
        File name of reference image.
    file_warp : str
        File name of initlal warp.
    file_jacobian : str
        Calculate and save jacobian of final warp field.
    file_out : str
        File name of output warp image.
    """
    # get filename
    path_out, _, _ = get_filename(file_out)

    # make output folder
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    command = "convertwarp"
    command += f" --ref={file_ref}"
    command += " --abs"
    command += f" --jacobian={file_jacobian}"
    command += " --relout"
    command += f" --warp1={file_warp}"
    command += f" --out={file_out}"

    # run
    execute_command(command)


def combine_warp(file_x, file_y, file_z, file_out):
    """Concatenate single warps in x-, y-, and z-direction into one warp file by
    contatenating in the fourth (time) dimension using FSL.

    Parameters
    ----------
    file_x : str
        File name of war in x-direction.
    file_y : str
        File name of warp in y-direction.
    file_z : str
        File name of warp in z-direction.
    file_out : str
        File name of merged warp file.
    """
    # get filename
    path_out, _, _ = get_filename(file_out)

    # make output folder
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    command = "fslmerge"
    command += " -t"
    command += f" {file_out}"
    command += f" {file_x} {file_y} {file_z}"

    # run
    execute_command(command)


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

    # run
    execute_command(command)
