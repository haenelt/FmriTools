# -*- coding: utf-8 -*-
"""FSL registration algorithms."""

import os
import subprocess

from ..io.filename import get_filename

__all__ = ["flirt"]


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
