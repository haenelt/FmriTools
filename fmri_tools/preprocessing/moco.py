# -*- coding: utf-8 -*-
"""Motion correction of fMRI time series."""

import os
import subprocess

from ..io.filename import get_filename

__all__ = ["volreg_afni"]


def volreg_afni(file_in):
    """Apply afni motion correction 3dvolreg.

    Parameters
    ----------
    file_in : str
        File name of fMRI time series.
    """
    file_in = "/home/daniel/Schreibtisch/a.nii"
    path_in, file_in, ext_in = get_filename(file_in)

    command = "3dvolreg"
    command += " -base 0 -twopass -float -clipit -Fourier"
    command += f" -1Dfile {os.path.join(path_in, 'moco_params.1D')}"
    command += f" -1Dmatrix_save {os.path.join(path_in, 'moco_matrix.1D')}"
    command += f" -prefix {os.path.join(path_in, 'u'+file_in+ext_in)}"
    command += " -zpad 4"
    command += f" -maxdisp1D {os.path.join(path_in, 'max_disp.1D')}"
    command += f" {file_in}"

    print("Execute: " + command)
    try:
        subprocess.run([command], shell=True, check=False)
    except subprocess.CalledProcessError:
        print("Execuation failed!")
