# -*- coding: utf-8 -*-
"""BOLD correction for VASO."""

import os
import subprocess

import nibabel as nb
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as Interp
from sh import gunzip

from ..io.filename import get_filename


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
