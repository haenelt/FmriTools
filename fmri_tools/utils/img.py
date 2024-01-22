# -*- coding: utf-8 -*-
"""Illustrate nifti images."""

import datetime
import os
import shutil as sh

import imageio as io
import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np

from ..segmentation.flat import orthographic_projection
from ..surface.filter import mris_fwhm

__all__ = ["get_gif", "get_movie", "get_retinotopy_images"]


def get_gif(img_file, path_output, name_output, nsteps, duration):
    """Create gif movie from n input images with fading from one image to the next
    image.

    Parameters
    ----------
    img_file : list
        List of input images.
    path_output : str
        Path where output is saved.
    name_output : str
        Base name of output file.
    nsteps : int
        Number of generated transition images.
    duration : float
        Duration for each frame in seconds.
    """
    # append first list item to the end of the list
    img_file.append(img_file[0])

    images = []
    for i in range(len(img_file) - 1):
        # input
        img1 = io.imread(img_file[i])
        img2 = io.imread(img_file[i + 1])

        # generates a list with transitions from img1 to img2
        for j in range(nsteps):
            img = np.multiply(img1, (nsteps - j) / nsteps) + np.multiply(
                img2, j / nsteps
            )
            img = img.astype("uint8")
            images.append(img)

    # save output gif
    io.mimwrite(
        os.path.join(path_output, name_output + ".gif"), images, duration=duration
    )


def get_movie(file_in, path_output, name_output, coord, axis=0, fps=10, transpose=True):
    """This function generates gif of a specific slice over a time series. It has the
    purpose to easily check preprocessing performance.

    Parameters
    ----------
    file_in : str
        File name of input file.
    path_output : str
        Path where output is saved.
    name_output : str
        Output file name without file extension.
    coord : int
        Selected slice.
    axis : int, optional
        Axis along which slice is selected. The default is 0.
    fps : float, optional
        Frames per second. The default is 10.
    transpose : bool, optional
        Rotate image. The default is True.

    Raises
    ------
    ValueError
        If `axis` number is invalid.
    """
    # make sub-folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # load time series
    data_img = nb.load(file_in)
    data_array = data_img.get_fdata()

    # get matrix dimension
    x_size = data_img.header["dim"][1]
    y_size = data_img.header["dim"][2]
    z_size = data_img.header["dim"][3]
    nt = data_img.header["dim"][4]

    # get movie array
    if axis == 0:
        movie_array = np.zeros([y_size, z_size, nt])
        for i in range(nt):
            movie_array[:, :, i] = data_array[coord, :, :, i]
    elif axis == 1:
        movie_array = np.zeros([x_size, z_size, nt])
        for i in range(nt):
            movie_array[:, :, i] = data_array[:, coord, :, i]
    elif axis == 2:
        movie_array = np.zeros([x_size, y_size, nt])
        for i in range(nt):
            movie_array[:, :, i] = data_array[:, :, coord, i]
    else:
        raise ValueError("Choose a valid axis!")

    # scale movie array
    movie_array = movie_array / np.max(movie_array) * 255
    movie_array = movie_array.astype("uint8")

    images = []
    for i in range(nt):
        images.append(movie_array[:, :, i])

    if transpose is True:
        images = np.transpose(images, axes=(0, 2, 1))

    # generate video file
    io.mimwrite(os.path.join(path_output, name_output + ".gif"), images, fps=fps)


def get_retinotopy_images(
    input_patch,
    input_vfs,
    input_phase,
    input_snr,
    input_white,
    hemi,
    path_output,
    img_res=0.2,
    theta=0,
    alpha=2,
    buffer=0,
    phase_fwhm=4,
    sigma=50,
    cleanup=True,
):
    """This function generates images for each step of a given retinotopy phase map from
    which animations of the temporal phase shift on the flattened surface can be made.

    Parameters
    ----------
    input_patch : str
        Filename of flattened patch.
    input_vfs : str
        Filename of visual fieldsign map.
    input_phase : str
        Filename of phase map.
    input_snr : str
        Filename of SNR map.
    input_white : str
        Filename of white surface.
    hemi : str
        Hemisphere.
    path_output : str
        Path where output is saved.
    img_res : float, optional
        Isotropic image resolution in mm. The default is 0.2.
    theta : float, optional
        Rotation of flat image in deg. The default is 0.
    alpha : float, optional
        Alpha shape value for concave hull computation. The default is 2.
    buffer : float, optional
        Smooth out concave hull.. The default is 0.
    phase_fwhm : float, optional
        Smoothing kernel for phase map smoothing in mm.. The default is 4.
    sigma : float, optional
        Gaussian kernel size for weighting towards single phase values.. The
        default is 50.
    cleanup : bool, optional
        Delete intermediate files. The default is True.

    Raises
    ------
    FileExistsError
        If the temporary folder already exists.
    """
    # phase step for single images
    phase_step = np.arange(-180, 181, 1)

    # set freesurfer path environment
    os.environ["SUBJECTS_DIR"] = path_output

    # freesurfer subject
    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = "".join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    sub = "tmp_" + tmp_string

    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # mimic freesurfer folder structure (with some additional folder for
    # intermediate files)
    path_sub = os.path.join(path_output, sub)
    path_img = os.path.join(path_output, "img")
    path_surf = os.path.join(path_sub, "surf")
    path_ortho = os.path.join(path_sub, "ortho")

    if not os.path.exists(path_sub):
        os.makedirs(path_sub)
    else:
        raise FileExistsError("Temporary folder already exists!")

    os.makedirs(path_surf)
    os.makedirs(path_ortho)
    os.makedirs(path_img)

    # copy surfaces to mimicked freesurfer folders
    sh.copyfile(input_patch, os.path.join(path_surf, hemi + ".patch"))
    sh.copyfile(input_white, os.path.join(path_surf, hemi + ".white"))
    sh.copyfile(input_vfs, os.path.join(path_surf, hemi + ".vfs.mgh"))
    sh.copyfile(input_phase, os.path.join(path_surf, hemi + ".phase.mgh"))
    sh.copyfile(input_snr, os.path.join(path_surf, hemi + ".snr.mgh"))

    # get orthographic projection
    orthographic_projection(
        os.path.join(path_surf, hemi + ".patch"),
        img_res,
        theta,
        alpha,
        buffer,
        path_ortho,
    )

    mris_fwhm(
        os.path.join(path_surf, hemi + ".phase.mgh"),
        os.path.join(path_surf, hemi + ".phase_smooth.mgh"),
        sub,
        phase_fwhm,
    )
    mris_fwhm(
        os.path.join(path_surf, hemi + ".snr.mgh"),
        os.path.join(path_surf, hemi + ".snr_smooth.mgh"),
        sub,
        phase_fwhm,
    )

    # read cmap, mask, vfs, phase
    cmap = (
        nb.load(os.path.join(path_ortho, hemi + ".patch.cmap.nii"))
        .get_fdata()
        .astype(int)
    )
    mask = nb.load(os.path.join(path_ortho, hemi + ".patch.mask.nii")).get_fdata()
    vfs = nb.load(os.path.join(path_surf, hemi + ".vfs.mgh")).get_fdata()
    phase = nb.load(os.path.join(path_surf, hemi + ".phase_smooth.mgh")).get_fdata()
    snr = nb.load(os.path.join(path_surf, hemi + ".snr_smooth.mgh")).get_fdata()

    # sample onto regular grid
    vfs_grid = np.zeros_like(cmap)
    phase_grid = np.zeros_like(cmap)
    snr_grid = np.zeros_like(cmap)
    for i in range(np.shape(cmap)[0]):
        for j in range(np.shape(cmap)[1]):
            if cmap[i, j] != 0:
                vfs_grid[i, j] = vfs[cmap[i, j]]
                phase_grid[i, j] = phase[cmap[i, j]]
                snr_grid[i, j] = snr[cmap[i, j]]

    # normalize snr
    snr_grid = snr_grid / np.max(snr_grid)

    # regular grid image of vfs
    img = np.zeros((np.shape(cmap)[0], np.shape(cmap)[1], 3))
    img[:, :, 0] = vfs_grid
    img[:, :, 1] = vfs_grid
    img[:, :, 2] = vfs_grid

    img[:, :, 0][img[:, :, 0] > -1] = 0
    img[:, :, 1][img[:, :, 1] > -1] = 0
    img[:, :, 2][img[:, :, 2] < 1] = 0
    img = np.abs(img)

    img[:, :, 0] = img[:, :, 0] * mask
    img[:, :, 1] = img[:, :, 1] * mask
    img[:, :, 2] = img[:, :, 2] * mask

    # regular grid images of single phases
    for k, _ in enumerate(phase_step):
        # get activation in red
        img2 = np.zeros((np.shape(cmap)[0], np.shape(cmap)[1], 4))

        img2[:, :, 0] = 1
        img2[:, :, 1] = 0
        img2[:, :, 2] = 0
        img2[:, :, 3] = _gaussian_filter(phase_grid, phase_step[k], sigma) * snr_grid

        img2[:, :, 0] = img2[:, :, 0] * mask
        img2[:, :, 1] = img2[:, :, 1] * mask
        img2[:, :, 2] = img2[:, :, 2] * mask

        # save plot
        fig = plt.figure(frameon=False)
        fig.set_size_inches(1, 1)
        ax = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
        ax.set_axis_off()
        fig.add_axes(ax)
        ax.imshow(img, aspect="auto")
        ax.imshow(img2, aspect="auto")
        fig.savefig(os.path.join(path_img, "img_" + str(k) + ".png"), dpi=400)
        plt.close("all")

    # delete intermediate files
    if cleanup:
        sh.rmtree(path_sub, ignore_errors=True)


def _gaussian_filter(x, x0, s):
    """Helper function for definition of gaussian filter.

    Parameters
    ----------
    x : (N, M) np.ndarray
        Input array.
    x0 : float
        Phase step.
    s : float
        Standard deviation of gaussian filter.

    Returns
    -------
    (N, M) np.ndarray
        Filtered array
    """
    g = 1 / (s * np.sqrt(2 * np.pi)) * np.exp(-((x - x0) ** 2) / (2 * s**2))
    return g / np.max(g)
