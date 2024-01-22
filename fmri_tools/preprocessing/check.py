# -*- coding: utf-8 -*-
"""Check preprocessing quality."""

import os
import shutil as sh
import subprocess

import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
from scipy.stats import pearsonr, shapiro

from ..io.filename import get_filename
from ..segmentation.mask import clean_ana, mask_ana, mask_epi
from ..utils.metrics import calc_mean

__all__ = ["check_preprocessing", "check_preprocessing_afni"]


# parameters for orig skullstrip
NITER_MASK = 3
SIGMA_MASK = 3


def check_preprocessing(files_in, file_mask, r_threshold=0.95, show_plots=False):
    """This scripts intends to identify corrupted volumes and runs in a set of
    functional time series to get a quantifiable and reproducible exclusion criterion
    (cf. Marquardt et al. 2017; Bergmann et al. 2019) by computing the spatial
    cross-correlation between signal intensities in functional time series. A regressor
    of no interest is written for each input time series to denote volumes below
    threshold.

    Parameters
    ----------
    files_in : list
        File name of preprocessed time series.
    file_mask : str
        File name of input mask.
    r_threshold : float, optional
        Pearson correlation threshold for outlier detection, by default 0.95.
    show_plots : bool, optional
        Show sanity plots, by default False.

    Returns
    -------
    None.

    """
    # get filename from first input entry
    _, name_file, _ = get_filename(files_in[0])

    # create output folder
    if len(files_in) < 2:
        path_output = os.path.join(os.path.dirname(files_in[0]), "correlation")
    else:
        path_output = os.path.join(
            os.path.dirname(os.path.dirname(files_in[0])), "correlation"
        )

    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # reference volumes is first volume of first time series
    data_0 = nb.load(files_in[0]).get_fdata()[:, :, :, 0]

    if file_mask is not None:
        mask_0 = nb.load(file_mask).get_fdata()

    # open logfile
    file = open(
        os.path.join(path_output, "correlation_" + name_file + ".txt"),
        "w",
        encoding="utf-8",
    )
    file.write("Percentage of volumes below threshold\n")
    file.write("Correlation threshold: " + str(r_threshold) + "\n\n")

    r_pearson_0 = []
    r_pearson = []
    r_shapiro = []
    p_pearson_0 = []
    p_pearson = []
    p_shapiro = []
    for i, file_ in enumerate(files_in):
        # load time series
        data_temp = nb.load(file_).get_fdata()

        # print progress
        print(f"Time series {i + 1}/{len(files_in)}")

        # create output folder for regressor of no interest
        path_logfile = os.path.join(os.path.dirname(file_), "outlier")
        if not os.path.exists(path_logfile):
            os.makedirs(path_logfile)

        # open logfile for regressor of no interest
        file2 = open(
            os.path.join(path_logfile, "correlation_regressor_" + name_file + ".txt"),
            "w",
            encoding="utf-8",
        )

        pearson_run = 0
        npearson_0 = 0
        nshapiro = 0
        for j in range(np.shape(data_temp)[3]):
            # load reference and current time step
            data_temp_0 = data_0.copy()
            data_temp_j = data_temp[:, :, :, j].copy()

            # mask time step
            if len(file_mask) > 0:
                data_temp_0 = data_temp_0[mask_0 == 1].flatten()
                data_temp_j = data_temp_j[mask_0 == 1].flatten()
            else:
                data_temp_0 = data_temp_0.flatten()
                data_temp_j = data_temp_j.flatten()

            # shapiro wilk
            r_tmp, p_tmp = shapiro(data_temp_j)
            r_shapiro = np.append(r_shapiro, r_tmp)
            p_shapiro = np.append(p_shapiro, p_tmp)

            if r_tmp < r_threshold:
                nshapiro += 1

            # pearson correlation to reference
            r_tmp, p_tmp = pearsonr(data_temp_0, data_temp_j)
            r_pearson_0 = np.append(r_pearson_0, r_tmp)
            p_pearson_0 = np.append(p_pearson_0, p_tmp)

            # sum for within-run correlation
            pearson_run += r_tmp

            if r_tmp < r_threshold:
                npearson_0 += 1
                file2.write("1\n")
            else:
                file2.write("0\n")

            if j < np.shape(data_temp)[3] - 1 and i < len(files_in):
                data_temp_1 = data_temp[:, :, :, j].copy()
                data_temp_2 = data_temp[:, :, :, j + 1].copy()

            elif j == np.shape(data_temp)[3] - 1 and i < len(files_in) - 1:
                data_temp_1 = data_temp[:, :, :, j].copy()
                data_temp_2 = nb.load(files_in[i + 1]).get_fdata()[:, :, :, 0]

            else:
                break

            # mask time step
            if len(file_mask) > 0:
                data_temp_1 = data_temp_1[mask_0 == 1].flatten()
                data_temp_2 = data_temp_2[mask_0 == 1].flatten()
            else:
                data_temp_1 = data_temp_1.flatten()
                data_temp_2 = data_temp_2.flatten()

            r_tmp, p_tmp = pearsonr(data_temp_1, data_temp_2)
            r_pearson = np.append(r_pearson, r_tmp)
            p_pearson = np.append(p_pearson, p_tmp)

        # close logfile for regressor of no interest
        file2.close()

        # percentage below threshold
        res_pearson_0 = npearson_0 / np.shape(data_temp)[3] * 100
        res_shapiro = nshapiro / np.shape(data_temp)[3] * 100

        # average within-run correlation
        pearson_run = pearson_run / np.shape(data_temp)[3]

        # update logfile
        file.write("Run: " + str(i + 1) + "\n")
        file.write("----------\n")
        file.write("Pearson (average within run): " + str(pearson_run) + "\n")
        file.write("Outlier percentage (pearson to ref): " + str(res_pearson_0) + "\n")
        file.write("Outlier percentage (shapiro): " + str(res_shapiro) + "\n\n\n")

    # close logfile
    file.close()

    # save variables
    np.savez(
        os.path.join(path_output, "correlation_" + name_file),
        r_shapiro=r_shapiro,
        p_shapiro=p_shapiro,
        r_pearson=r_pearson,
        p_pearson=p_pearson,
        r_pearson_0=r_pearson_0,
        p_pearson_0=p_pearson_0,
    )

    # plots
    fig, ax = plt.subplots()
    ax.plot(r_shapiro, "r")
    ax.set_xlabel("Time in TR")
    ax.set_ylabel("r-value (Shapiro-Wilk test)")
    ax.set_title("Shapiro-Wilk test for each time step")
    ax.hlines(r_threshold, 0, len(r_shapiro), linestyle="dashed")
    fig.savefig(
        os.path.join(path_output, "r_shapiro_" + name_file + ".png"),
        format="png",
        bbox_inches="tight",
    )
    if show_plots:
        plt.show()

    fig, ax = plt.subplots()
    ax.plot(p_shapiro, "r")
    ax.set_xlabel("Time in TR")
    ax.set_ylabel("p-value (Shapiro-Wilk test)")
    ax.set_title("Shapiro-Wilk test for each time step")
    fig.savefig(
        os.path.join(path_output, "p_shapiro_" + name_file + ".png"),
        format="png",
        bbox_inches="tight",
    )
    if show_plots:
        plt.show()

    fig, ax = plt.subplots()
    ax.plot(r_pearson_0, "r")
    ax.set_xlabel("Time in TR")
    ax.set_ylabel("r-value (Pearson correlation)")
    ax.set_title("Time series spatial correlation to volume ref")
    ax.hlines(r_threshold, 0, len(r_pearson_0), linestyle="dashed")
    fig.savefig(
        os.path.join(path_output, "r_pearson_0_" + name_file + ".png"),
        format="png",
        bbox_inches="tight",
    )
    if show_plots:
        plt.show()

    fig, ax = plt.subplots()
    ax.plot(p_pearson_0, "r")
    ax.set_xlabel("Time in TR")
    ax.set_ylabel("p-value (Pearson correlation)")
    ax.set_title("Time series spatial correlation to volume ref")
    fig.savefig(
        os.path.join(path_output, "p_pearson_0_" + name_file + ".png"),
        format="png",
        bbox_inches="tight",
    )
    if show_plots:
        plt.show()

    fig, ax = plt.subplots()
    ax.plot(r_pearson, "r")
    ax.set_xlabel("Time in TR")
    ax.set_ylabel("r-value (Pearson correlation)")
    ax.set_title("Time series spatial correlation to volume i-1")
    ax.hlines(r_threshold, 0, len(r_pearson), linestyle="dashed")
    fig.savefig(
        os.path.join(path_output, "r_pearson_" + name_file + ".png"),
        format="png",
        bbox_inches="tight",
    )
    if show_plots:
        plt.show()

    fig, ax = plt.subplots()
    ax.plot(p_pearson, "r")
    ax.set_xlabel("Time in TR")
    ax.set_ylabel("p-value (Pearson correlation)")
    ax.set_title("Time series spatial correlation to volume i-1")
    fig.savefig(
        os.path.join(path_output, "p_pearson_" + name_file + ".png"),
        format="png",
        bbox_inches="tight",
    )
    if show_plots:
        plt.show()


def check_preprocessing_afni(
    file_epi, file_t1, file_mask, qthr=0.001, pthr=0.01, cleanup=False
):
    """This script looks for outliers in a functional time series using the Afni
    function 3dToutcount. The fraction of outliers found within a volume defined by a
    mask are written out for each time point. Additionally, a graphical visualization, a
    regressor of no interest and a short summary are written. The script needs an
    installation of fsl and afni.

        Parameters
        ----------
        file_epi : str
            File name of preprocessed time series.
        file_t1 : str
            File name of anatomical image.
        file_mask : str
            File name of mask in space of anatomical image.
        qthr : float, optional
            Threshold defined by AFNI, by default 0.001.
        pthr : float, optional
            Outlier threshold, by default 0.01.
        cleanup : bool, optional
            Remove intermediate files, by default False.
    """
    # get fileparts from first input entry
    path_epi, name_epi, ext_epi = get_filename(file_epi)
    _, _, ext_t1 = get_filename(file_t1)
    _, _, ext_mask = get_filename(file_mask)

    # make folder structure
    path_output = os.path.join(path_epi, "outlier_afni")
    path_temp = os.path.join(path_output, "temp")

    if not os.path.exists(path_output):
        os.makedirs(path_output)

    if not os.path.exists(path_temp):
        os.makedirs(path_temp)

    # copy input files into temporary folder
    sh.copyfile(file_t1, os.path.join(path_temp, "T1" + ext_t1))
    sh.copyfile(file_mask, os.path.join(path_temp, "mask" + ext_mask))

    # get mean time series
    calc_mean(file_epi, path_temp, "epi", method="mean")

    # new filenames
    file_epi_mean = os.path.join(path_temp, "mean_epi" + ext_epi)
    file_t1_copy = os.path.join(path_temp, "T1" + ext_t1)
    file_mask_copy = os.path.join(path_temp, "mask" + ext_mask)

    # get mask
    clean_ana(file_t1_copy, 1000.0, 4095.0, overwrite=True)  # clean ana
    mask_ana(file_t1_copy, file_mask_copy, background_bright=False)  # mask t1
    mask_epi(
        file_epi_mean,
        os.path.join(path_temp, "pT1" + ext_t1),
        file_mask_copy,
        NITER_MASK,
        SIGMA_MASK,
    )

    # get outlier count within mask
    command = "3dToutcount"
    command += f" -mask {os.path.join(path_temp, 'mask_def2.nii.gz')}"
    command += " -fraction"
    command += f" -qthr {qthr}"
    command += f" {file_epi}"
    command += f" > {os.path.join(path_output, 'outlier_afni.txt')}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")

    # make plot
    log_data = np.loadtxt(os.path.join(path_output, "outlier_afni.txt"))

    # plots
    fig, ax = plt.subplots()
    ax.plot(log_data, "r")
    ax.set_xlabel("Time in TR")
    ax.set_ylabel("Outlier count")
    ax.set_title("Fraction of outliers found at each time step")
    ax.hlines(pthr, 0, len(log_data), linestyle="dashed")
    fig.savefig(
        os.path.join(path_output, "outlier_afni.png"), format="png", bbox_inches="tight"
    )

    # make regressor of no interest
    log_outlier = np.zeros_like(log_data)
    log_outlier[log_data > pthr] = 1
    np.savetxt(
        os.path.join(path_output, "outlier_afni_regressor_" + name_epi + ".txt"),
        log_outlier,
        fmt="%d",
    )

    # make summary textfile
    file = open(
        os.path.join(path_output, "outlier_afni_summary.txt"), "w", encoding="utf-8"
    )
    file.write("qthr: " + str(qthr) + "\n")
    file.write("pthr: " + str(pthr) + "\n")
    file.write("Number of outliers: " + str(np.sum(log_outlier)))
    file.close()

    # clean intermediate files
    if cleanup:
        sh.rmtree(path_temp, ignore_errors=True)
