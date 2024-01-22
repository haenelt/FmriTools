# -*- coding: utf-8 -*-
"""Surface smoothing."""

import os
import subprocess

from ..io.filename import get_filename

__all__ = ["mris_smooth"]


def mris_smooth(file_in, file_out, n_iter):
    """Smooth surface tesselation.

    Parameters
    ----------
    file_in : str
        File name of input surface.
    file_out : str
        File name of output surface,
    n_iter : int
        Number of smoothing iterations.
    """
    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    command = "mris_smooth -nw"
    command += f" -n {n_iter}"
    command += f" {file_in}"
    command += f" {file_out}"

    print("Execute: " + command)
    try:
        subprocess.run([command], shell=True, check=False)
    except subprocess.CalledProcessError:
        print("Execuation failed!")
