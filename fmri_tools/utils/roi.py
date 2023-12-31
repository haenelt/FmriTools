# -*- coding: utf-8 -*-
"""Region of interest (ROI) manipulations."""

import os
import subprocess

from ..io.filename import get_filename

__all__ = ["extract_vol"]


def extract_vol(file_in, file_out, t_min, t_size):
    """Extract 3D volumes from nifti time series.

    Parameters
    ----------
    file_in : str
        File name of image time series.
    file_out : str
        File name of extracted time series.
    t_min : int
        Start volume.
    t_size : int
        Number of volumes.
    """
    # get filename
    _, _, ext_in = get_filename(file_in)
    path_out, _, ext_out = get_filename(file_out)

    if ext_in or ext_out not in [".nii", ".nii.gz"]:
        raise ValueError("Invalid file extension!")

    # make output folder
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    command = "fslroi"
    command += f" {file_in}"
    command += f" {file_out}"
    command += f" {t_min}"
    command += f" {t_size}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")
