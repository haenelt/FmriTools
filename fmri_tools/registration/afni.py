# -*- coding: utf-8 -*-
"""AFNI realignment commands."""

from pathlib import Path
from fmri_tools import execute_command
from fmri_tools.io.vol import gzip_nii, gunzip_nii

__all__ = ["volreg", "allineate", "prepare_header", "extract_ref"]


def prepare_header(file_in):
    """This function handles potential oblique data sets, which can affect motion
    correction between across datasets, e.g. run, if not handled beforehand. Be aware
    that this function overwrites the input dataset!!!

    Since nifti_tool cannot process compressed datasets, compressed input files are
    temporarily uncompressed before processing.

    Parameters
    ----------
    file_in : str
        File name of input time series.
    """
    ext = "".join(Path(file_in).suffixes)
    uncompress = False
    if ".gz" in ext:
        uncompress = True
        file_in = gunzip_nii(file_in, overwrite=True, remove_original=True)

    # overwrite dataset as afni saves data
    command = "3dcalc"
    command += f" -a {file_in}"
    command += " -expr 'a'"
    command += f" -prefix {file_in}"
    command += " -overwrite"

    # run
    execute_command(command)

    # handle affine matrices
    command = "nifti_tool"
    command += " -mod_hdr -mod_field sform_code 1 -mod_field qform_code 1"
    command += f" -infiles {file_in}"
    command += " -overwrite"

    # run
    execute_command(command)

    # deoblique input data set
    command = "3drefit"
    command += " -deoblique"
    command += f" {file_in}"

    # run
    execute_command(command)

    if uncompress:
        gzip_nii(file_in, overwrite=True, remove_original=True)


def extract_ref(file_in, file_out, t=4):
    """Extract a volume from an fMRI time series, which can be used as reference volume
    in motion correction usign afni.

    file_in : str
        File name of input time series.
    file_out : str
        File name of reference volume.
    t : int
        Which time point in the time series to extract. Per default, not the first
        volume is selected to be sure that steady-state magnetization is reached.
    """
    dir_out = Path(file_out).parent
    dir_out.mkdir(exist_ok=True, parents=True)

    command = "3dbucket"
    command += f" -prefix {file_out}"
    command += f" {file_in}'[{t}]'"

    # run
    execute_command(command)


def volreg(file_in, file_out, file_ref):
    """Apply afni motion correction 3dvolreg.

    Parameters
    ----------
    file_in : str
        File name of input time series.
    file_out : str
        File name of realigned output time series.
    file_ref : str
        File name of target volume.
    """
    dir_out = Path(file_out).parent
    dir_out.mkdir(exist_ok=True, parents=True)

    command = "3dvolreg"
    command += f" -base {file_ref}"
    command += " -twopass -float -clipit -Fourier"
    command += f" -1Dfile {dir_out / 'moco_params.1D'}"
    command += f" -1Dmatrix_save {dir_out / 'moco_matrix.1D'}"
    command += f" -prefix {file_out}"
    command += " -zpad 4"
    command += f" -maxdisp1D {dir_out / 'max_disp.1D'}"
    command += f" {file_in}"

    # run
    execute_command(command)


def allineate(file_in, file_out, file_ref, file_moco):
    """Apply motion estimates to input data using afni.

    Parameters
    ----------
    file_in : str
        File name of input time series.
    file_out : str
        File name of realigned output time series.
    file_ref : str
        File name of target volume.
    file_moco : str
        File name of realignment motion paraterms as saved by afni's 3dvolreg in *.1D
        format.
    """
    dir_out = Path(file_out).parent
    dir_out.mkdir(exist_ok=True, parents=True)

    command = "3dAllineate"
    command += f" -prefix {file_out}"
    command += f" -1Dmatrix_apply {file_moco}"
    command += " -final wsinc5"
    command += " -cubic -float"
    command += f" -base {file_ref}"
    command += f" -input {file_in}"

    # run
    execute_command(command)
