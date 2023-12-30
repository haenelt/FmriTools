# -*- coding: utf-8 -*-

import datetime
import os
import shutil as sh

import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np

from ..segmentation.orthographic_projection import orthographic_projection


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
    """Get retinotopy images.

    This function generates images for each step of a given retinotopy phase map
    from which animations of the temporal phase shift on the flattened surface
    can be made.

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

    Returns
    -------
    None.

    """

    def gaussian_filter(x, x0, s):
        # return gaussian filter
        g = 1 / (s * np.sqrt(2 * np.pi)) * np.exp(-((x - x0) ** 2) / (2 * s**2))
        g = g / np.max(g)

        return g

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

    os.system(
        "mris_fwhm"
        + " --s "
        + sub
        + " --hemi "
        + hemi
        + " --smooth-only "
        + " --fwhm "
        + str(phase_fwhm)
        + " --i "
        + os.path.join(path_surf, hemi + ".phase.mgh")
        + " --o "
        + os.path.join(path_surf, hemi + ".phase_smooth.mgh")
    )

    os.system(
        "mris_fwhm"
        + " --s "
        + sub
        + " --hemi "
        + hemi
        + " --smooth-only "
        + " --fwhm "
        + str(phase_fwhm)
        + " --i "
        + os.path.join(path_surf, hemi + ".snr.mgh")
        + " --o "
        + os.path.join(path_surf, hemi + ".snr_smooth.mgh")
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
    for k in range(len(phase_step)):
        # get activation in red
        img2 = np.zeros((np.shape(cmap)[0], np.shape(cmap)[1], 4))

        img2[:, :, 0] = 1
        img2[:, :, 1] = 0
        img2[:, :, 2] = 0
        img2[:, :, 3] = gaussian_filter(phase_grid, phase_step[k], sigma) * snr_grid

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
