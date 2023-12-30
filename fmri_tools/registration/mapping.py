# -*- coding: utf-8 -*-
"""Map to surface."""

import subprocess

from ..io.filename import get_filename

__all__ = ["mri_vol2surf"]


def mri_vol2surf(file_in, file_out, sub, interp_method="nearest"):
    """Use freesurfer mri_vol2surf to sample volume data onto a surface mesh.

    Parameters
    ----------
    file_in : str
        File name of input volume file.
    file_out : str
        File name of output overlay file.
    sub : str
        Name of freesurfer subject.
    interp_method : str, optional
        Interpolation method (nearest or trilinear), by default "nearest"
    """
    # get hemisphere
    _, _, hemi = get_filename(file_in)
    hemi = hemi.replace(".", "")
    if hemi not in ["lh", "rh"]:
        raise ValueError("No hemisphere specified in file name!")

    command = "mri_vol2surf"
    command += f" --hemi {hemi}"
    command += f" --interp {interp_method}"
    command += f" --o {file_out}"
    command += " --out_type mgh"
    command += f" --regheader {sub}"
    command += " --projdist 0.000"
    command += f" --mov {file_in}"
    command += " --surf source"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")
