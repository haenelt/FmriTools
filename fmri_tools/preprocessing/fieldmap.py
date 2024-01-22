# -*- coding: utf-8 -*-
"""Surface smoothing."""

import os
import subprocess

from ..io.filename import get_filename

__all__ = ["prepare_fieldmap_fsl", "fugue_fsl"]


def prepare_fieldmap_fsl(file_magn, file_phase, file_out, delta_te):
    """Create fieldmap from siemens phase difference image using FSL.

    Parameters
    ----------
    file_magn : str
        File name of magnitude image (take first echo).
    file_phase : str
        File name of phase differene image.
    file_out : str
        File name of fieldmap image (in rad/s).
    delta_te : float
        Delta TE in ms.

    Returns
    -------
    None.

    """
    # make output folder
    path_out, _, _ = get_filename(file_out)
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    command = "fsl_prepare_fieldmap SIEMENS"
    command += f" {file_phase}"
    command += f" {file_magn}"
    command += f" {file_out}"
    command += f" {delta_te:.6f}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")


def fugue_fsl(file_in, file_fmap, file_shift, dwell_time, gaussian_sigma, udir):
    """Unwarp image based on field map. The unwarped input file is saved in the same
    directory.

    Parameters
    ----------
    file_in : str
        File name of input image.
    file_fmap : str
        File name of field map (in rad/s).
    file_shift : str
        File name of output image containing applied shift.
    dwell_time : float
        Dwell time (in s).
    gaussian_sigma : float
        Apply Gaussian smoothing of sigma (in mm).
    udir : str
        Unwarping direction (x, y, z, x-, y-, or z-).

    Returns
    -------
    None.

    """
    # output file name
    _, name_in, ext_in = get_filename(file_in)

    # make output folder
    path_out, _, _ = get_filename(file_shift)
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    command = "fugue"
    command += f" --dwell={dwell_time:.10f}"
    command += f" --loadfmap={file_fmap}"
    command += f" --in={file_in}"
    command += f" --saveshift={file_shift}"
    command += f" --smooth3={gaussian_sigma:.2f}"
    command += f" --unwarpdir={udir}"
    command += f" --unwarp={name_in}_unwarped{ext_in}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")
