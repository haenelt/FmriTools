# -*- coding: utf-8 -*-
"""BOLD correction for VASO."""

import os
import subprocess

import nibabel as nb
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as Interp
from sh import gunzip

from ..io.filename import get_filename

__all__ = [
    "split_vaso",
    "boco_vaso_regrid",
    "boco_vaso",
    "regrid_time_series",
    "regrid_time_series_afni",
]


def split_vaso(file_in, start_vol=2, end_vol=5, first_bold=True):
    """This function splits the vaso and bold time series into separate files. Splitted
    time series are names as bold (not-nulled) and vaso (nulled), respectively.
    Additionally, volumes at the beginning (non steady-state volumes) can be overwritten
    by following steady-state volumes. Furthermore, volumes at the end of the scan can
    be discarded. Note that the number of volumes refers here to the individual volumes
    in the bold+vaso time series.

    Parameters
    ----------
    file_in : _type_
        File name of SS-SI VASO time series.
    start_vol : int, optional
        Overwrite volumes at the beginning, by default 2.
    end_vol : int, optional
        Discard volumes at the end, by default 5.
    first_bold : bool, optional
        If True, time series starts with (not-nulled) bold volume, by default True.

    Returns
    -------
    None.

    """
    # load data
    data = nb.load(file_in)
    data_array = data.get_fdata()

    # overwrite non steady-state volumes
    data_array[:, :, :, 0:start_vol] = data_array[:, :, :, start_vol : 2 * start_vol]

    # discard volumes at the end
    if end_vol != 0:
        data_array = data_array[:, :, :, :-end_vol]

    # split into even and odd runs
    t_even = np.arange(0, np.shape(data_array)[3], 2)
    t_odd = np.arange(1, np.shape(data_array)[3], 2)

    if first_bold:
        bold_array = data_array[:, :, :, t_even]
        vaso_array = data_array[:, :, :, t_odd]
    else:
        bold_array = data_array[:, :, :, t_odd]
        vaso_array = data_array[:, :, :, t_even]

    # new array length
    data.header["dim"][4] = np.shape(vaso_array)[3]

    output = nb.Nifti1Image(vaso_array, data.affine, data.header)
    nb.save(output, os.path.join(os.path.dirname(file_in), "vaso.nii"))

    output = nb.Nifti1Image(bold_array, data.affine, data.header)
    nb.save(output, os.path.join(os.path.dirname(file_in), "bold.nii"))


def boco_vaso(file_vaso, file_bold, TR, vaso_threshold=6.0):
    """This scripts corrects a vaso time series for bold contamination. First, both time
    series are upsampled and the vaso time series is shifted by one time step. BOLD
    correction is performed by dividing both time series. In the end, unrealistic vaso
    values are removed. The script needs an installation of afni.

    Parameters
    ----------
    file_vaso : str
        File name of nulled time series.
    file_bold : str
        File name of not-nulled time series.
    TR : float
        Repetition time of single volumes in s.
    vaso_threshold : float, optional
        Threshold unrealistic values, by default 6.

    Returns
    -------
    None.

    """
    # prepare path and filename
    path_vaso, name_vaso, ext_vaso = get_filename(file_vaso)
    path_bold, name_bold, ext_bold = get_filename(file_bold)

    file_vaso_upsampled = os.path.join(path_vaso, f"{name_vaso}_upsampled.{ext_vaso}")
    file_bold_upsampled = os.path.join(path_bold, f"{name_bold}_upsampled.{ext_bold}")

    # upsample vaso and bold time series
    command = "3dUpsample -overwrite -datum short"
    command += f" -prefix {file_vaso_upsampled}"
    command += " -n 2"
    command += f" -input {file_vaso}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")

    command = "3dUpsample -overwrite -datum short"
    command += f" -prefix {file_bold_upsampled}"
    command += " -n 2"
    command += f" -input {file_bold}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")

    # load vaso data and shift in time
    vaso = nb.load(file_vaso_upsampled)
    vaso_array = vaso.get_fdata()
    vaso_array = vaso_array[:, :, :, :-1]
    vaso_array = np.concatenate(
        (np.expand_dims(vaso_array[:, :, :, 0], axis=3), vaso_array), axis=3
    )

    # load bold data
    bold = nb.load(file_bold_upsampled)
    bold_array = bold.get_fdata()

    # bold correction
    vaso_array = np.divide(vaso_array, bold_array)

    # remove nans and infs
    vaso_array[np.isnan(vaso_array)] = 0
    vaso_array[np.isinf(vaso_array)] = 0

    # clean vaso data that are unrealistic
    vaso_array[vaso_array < 0] = 0
    vaso_array[vaso_array >= vaso_threshold] = vaso_threshold

    output = nb.Nifti1Image(vaso_array, vaso.affine, vaso.header)
    file_vaso_corrected = os.path.join(
        path_vaso, f"{name_vaso}_upsampled_corrected.{ext_vaso}"
    )
    nb.save(output, file_vaso_corrected)

    # change TR in header
    command = "3drefit"
    command += f" -TR {TR}"
    command += f" {file_bold_upsampled}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")

    command = "3drefit"
    command += f" -TR {TR}"
    command += f" {file_vaso_upsampled}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")

    command = "3drefit"
    command += f" -TR {TR}"
    command += f" {file_vaso_corrected}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")


def boco_vaso_regrid(
    file_vaso,
    file_bold,
    TR_old,
    TR_new,
    start_bold,
    start_vaso,
    nvol_remove=0,
    vaso_threshold=6.0,
):
    """This scripts corrects a vaso time series for bold contamination. First, both time
    series are upsampled to a common time grid. BOLD correction is performed by dividing
    both time series. In the end, unrealistic vaso values are removed. The script needs
    an installation of afni.

    Parameters
    ----------
    file_vaso : str
        File name of nulled time series.
    file_bold : str
        File name of not-nulled time series.
    TR_old : float
        Original effective repetition time in s (nulled + not-nulled).
    TR_new : float
        New repetition time after BOLD correction in s.
    start_bold : float
        Start of not-nulled block in s.
    start_vaso : float
        Start of nulled block in s.
    nvol_remove : int, optional
        Number of volumes removed at the end of the time series, by default 0.
    vaso_threshold : float, optional
        Threshold unrealistic values, by default 6.

    Returns
    -------
    None.

    """
    # get filenames
    path_bold, name_bold, ext_bold = get_filename(file_bold)
    path_vaso, name_vaso, ext_vaso = get_filename(file_vaso)

    # upsample time series
    regrid_time_series(
        file_bold,
        path_bold,
        TR_old,
        TR_new,
        t_start=start_bold,
        nvol_remove=nvol_remove,
    )
    regrid_time_series(
        file_vaso,
        path_vaso,
        TR_old,
        TR_new,
        t_start=start_vaso,
        nvol_remove=nvol_remove,
    )

    # new filenames
    file_bold = os.path.join(path_bold, name_bold + "_upsampled" + ext_bold)
    file_vaso = os.path.join(path_vaso, name_vaso + "_upsampled" + ext_vaso)
    file_vaso_corrected = os.path.join(
        path_vaso, name_vaso + "_upsampled_corrected" + ext_vaso
    )

    # load bold data
    bold_data = nb.load(file_bold)
    bold_array = bold_data.get_fdata()

    # load vaso data
    vaso_data = nb.load(file_vaso)
    vaso_array = vaso_data.get_fdata()

    # bold correction
    vaso_array = np.divide(vaso_array, bold_array)

    # remove nans and infs
    vaso_array[np.isnan(vaso_array)] = 0
    vaso_array[np.isinf(vaso_array)] = 0

    # clean vaso data that are unrealistic
    vaso_array[vaso_array < 0] = 0
    vaso_array[vaso_array >= vaso_threshold] = vaso_threshold

    # write output
    output = nb.Nifti1Image(vaso_array, vaso_data.affine, vaso_data.header)
    nb.save(output, file_vaso_corrected)


def regrid_time_series(file_in, path_output, tr_old, tr_new, t_start=0, nvol_remove=0):
    """This function interpolates the time series onto a new time grid using cubic
    interpolation. Only for writing the new TR in the header of the output time series,
    AFNI has to be included in the search path.

    Parameters
    ----------
    file_in : str
        Time series filename.
    path_output : str
        Path where output is written.
    tr_old : float
        TR of time series in s.
    tr_new : float
        TR of regridded time series in s.
    t_start : float, optional
        Shift time series in s (t_start >= 0 and <= TR_old). The default is 0.
    nvol_remove : int, optional
        Remove volumes at the end of the time series.

    Returns
    -------
    None.

    """
    # get filename
    _, name_input, ext_input = get_filename(file_in)

    # print to console
    print("time series regridding for: " + name_input)

    # load data
    data = nb.load(file_in)
    nx = data.header["dim"][1]
    ny = data.header["dim"][2]
    nz = data.header["dim"][3]
    nt = data.header["dim"][4]

    # get time grid
    tt = tr_old * nt  # total acquisition time
    tr_append = (
        np.floor(t_start / tr_old + 1).astype(int) * tr_old
    )  # number of appended TRs in input array

    # input grid
    t_old = np.arange(-tr_append, tt + tr_append, tr_old) + t_start

    # output grid
    t_new = np.arange(0, tt + tr_append, tr_new)
    t_new_append = np.flip(np.arange(0, -tr_append, -tr_new)[1:])
    if not len(t_new_append):
        t_new_append = -tr_new
    t_new = np.append(t_new_append, t_new)

    # load array with appended volumes
    n_append = int(tr_append / tr_old)
    data_array = np.zeros((nx, ny, nz, nt + 2 * n_append))
    data_array[:, :, :, n_append:-n_append] = data.get_fdata()
    for i in range(n_append):
        data_array[:, :, :, i] = data.get_fdata()[:, :, :, 0]
        data_array[:, :, :, -(i + 1)] = data.get_fdata()[:, :, :, -1]

    # temporal interpolation
    data_array_regrid = np.zeros((nx, ny, nz, len(t_new)))
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                cubic_interper = Interp(t_old, data_array[x, y, z, :], k=3)
                data_array_regrid[x, y, z, :] = cubic_interper(t_new)

    # delete appended volumes
    vols_keep1 = t_new >= 0
    vols_keep2 = t_new < tt
    vols_keep = vols_keep1 * vols_keep2
    data_array_regrid = data_array_regrid[:, :, :, vols_keep]

    # clean corrected array
    data_array_regrid[np.isnan(data_array_regrid)] = 0
    data_array_regrid[data_array_regrid < 0] = 0

    # remove volumes at the end
    if nvol_remove:
        data_array_regrid = data_array_regrid[:, :, :, :-nvol_remove]

    # update data header
    data.header["dim"][4] = np.shape(data_array_regrid)[3]
    data.header["datatype"] = 16

    # write output
    output = nb.Nifti1Image(data_array_regrid, data.affine, data.header)
    nb.save(output, os.path.join(path_output, name_input + "_upsampled" + ext_input))

    # change TR in header
    command = "3drefit"
    command += f" -TR {tr_new}"
    command += f" {os.path.join(path_output, name_input + '_upsampled' + ext_input)}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")


def regrid_time_series_afni(file_in, n=2):
    """This function upsamples a time series using the afni function 3dUpsample. Before
    running the function, set the afni environment by calling AFNI in the terminal.
    Output is an upsampled nifti time series.

    Parameters
    ----------
    file_in : str
        Time series filename.
    n : int, optional
        Upsampling factor. The default is 2.

    Returns
    -------
    None.

    """
    clean_unzip = 0
    if os.path.splitext(file_in)[1] == ".gz":
        gunzip(file_in)
        clean_unzip = 1
        file_in = os.path.splitext(file_in)[0]

    # prepare path and filename
    path_file = os.path.dirname(file_in)
    name_file = os.path.splitext(os.path.basename(file_in))[0]

    # upsample vaso and bold time series
    command = "3dUpsample -overwrite -datum short "
    command += f" -prefix {os.path.join(path_file, name_file + '_upsampled.nii')}"
    command += f" -n {n}"
    command += f" -input {file_in}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")

    if clean_unzip:
        os.remove(file_in)
