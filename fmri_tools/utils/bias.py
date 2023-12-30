# -*- coding: utf-8 -*-
"""Handling spatial bias fields."""

import os
import subprocess

import nibabel as nb
import numpy as np

from ..io.filename import get_filename

__all__ = ["remove_bias_ants", "robust_combination"]


def remove_bias_ants(file_in, file_out, save_bias=True):
    """Apply the ANTS N4 bias field correction.

    Parameters
    ----------
    file_in : str
        File name of input image.
    file_out : str
        File name of output image.
    save_bias : bool, optional
        Save bias field to disk.
    """
    # file name of saved bias field
    path, _, _ = get_filename(file_out)
    file_out_bias = os.path.join(path, "n4bias.nii")

    command = "N4BiasFieldCorrection"
    command += " -d 3"
    command += f" --input-image {file_in}"
    if save_bias:
        command += f" --output [ {file_out}, {file_out_bias}]"
    else:
        command += f" --output {file_out}"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")


def robust_combination(fileUNI, fileINV1, fileINV2, regularisation, path_output):
    """This script implements the regularisation proposed in [1], which allows the
    creation of MP2RAGE flat images without the strong background noise in air regions.

    Although in the original paper the method only worked on raw multichannel data, here
    that constrain has been overcome and the correction can be implemented if both SOS
    images of the two inversion times exist and a MP2RAGE flat image that has been
    calculated directly from the multichannel data as initially proposed in [2].

    This script is a tranlation from Marques Matlab script found on his github
    repository [3].

    The output should ideally look like the standard MPRAGE (no noise in the
    background). If it has too much noise on the background, give it a bigger
    regularisation value.

    If the value is too big you will start noticing that the image gets a bigger bias
    field - which will deteriorate the segmentation results. Once you are happy with the
    value you found for one subject you can use the same for all the following subjects
    by just calling the function with the same regularisation parameter.

    Usually the regularisation factor shouldn't be greater than 10, but that is not the
    case when the image is bias field corrected, in which case the noise estimates at
    the edge of the image <noiselevel> might not be such a good measure.

    Note: To omit the warnings by running the rootsquares help functions, I set all
    voxel values with 0 to nan. After running that function, these voxels are set back
    to 0.

    Parameters
    ----------
    fileUNI : str
        Path of flat image.
    fileINV1 : str
        Path of corrsponding first inversion recovery.
    fileINV2 : str
        Path of corresponding second inversion recovery.
    regularisation : float
        Regularisation parameter for background removal.
    path_output : str
        Path where all output images are saved.

    Returns
    -------
    None.

    References
    -------
    .. [1] O'Brien, KR, et al, Robust T1-weighted structural brain imaging and
    morphometry at 7T using MP2RAGE, PLoS ONE 9(6), 1--7 (2014).
    .. [2] Marques, JP, et al. MP2RAGE, a self bias-field corrected sequence for
    improved segmentation and T1-mapping at high field, Neuroimage 49(2),
    1271--1281 (2010).
    .. [3] https://github.com/JosePMarques/MP2RAGE-related-scripts
    (accessed 02-10-2018)

    """
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # define relevant functions
    def MP2RAGErobustfunc(INV1, INV2, beta):
        return (np.conj(INV1) * INV2 - beta) / (INV1**2 + INV2**2 + 2 * beta)

    def rootsquares_pos(a, b, c):
        return (-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)

    def rootsquares_neg(a, b, c):
        return (-b - np.sqrt(b**2 - 4 * a * c)) / (2 * a)

    # load data
    MP2RAGE_img = nb.load(fileUNI)
    MP2RAGE_array = MP2RAGE_img.get_fdata()
    INV1_img = nb.load(fileINV1)
    INV1_array = INV1_img.get_fdata()
    INV2_img = nb.load(fileINV2)
    INV2_array = INV2_img.get_fdata()

    # convert MP2RAGE to [-0.5,0.5] if not in this range (assumes that only
    # positive values are found)
    if np.min(MP2RAGE_array) >= 0 and np.max(MP2RAGE_array) >= 0.51:
        MP2RAGE_array = (MP2RAGE_array - np.max(MP2RAGE_array) / 2) / np.max(
            MP2RAGE_array
        )
        integerformat = 1
    else:
        integerformat = 0

    # MP2RAGE is a phase sensitive coil combination.. some more maths has to be
    # performed to get a better INV1 estimate which here is done by assuming
    # both INV2 is closer to a real phase sensitive combination

    # give the correct polarity to INV1 (phase correction)
    INV1_array = np.sign(MP2RAGE_array) * INV1_array

    # change all values in MP2RAGE_array = 0 to nan
    MP2RAGE_array[MP2RAGE_array == 0] = np.nan

    INV1pos_array = rootsquares_pos(
        -MP2RAGE_array, INV2_array, -(INV2_array**2) * MP2RAGE_array
    )
    INV1neg_array = rootsquares_neg(
        -MP2RAGE_array, INV2_array, -(INV2_array**2) * MP2RAGE_array
    )

    # change all nan back to 0
    INV1_array[np.isnan(INV1_array)] = 0
    INV1pos_array[np.isnan(INV1pos_array)] = 0
    INV1neg_array[np.isnan(INV1neg_array)] = 0

    INV1final_array = INV1_array.copy()
    INV1final_array[
        np.abs(INV1_array - INV1pos_array) > np.abs(INV1_array - INV1neg_array)
    ] = INV1neg_array[
        np.abs(INV1_array - INV1pos_array) > np.abs(INV1_array - INV1neg_array)
    ]
    INV1final_array[
        np.abs(INV1_array - INV1pos_array) <= np.abs(INV1_array - INV1neg_array)
    ] = INV1pos_array[
        np.abs(INV1_array - INV1pos_array) <= np.abs(INV1_array - INV1neg_array)
    ]

    # lambda calculation
    noiselevel = regularisation * np.mean(INV2_array[:, -10:, -10:])
    MP2RAGERobustPhaseSensitive_array = MP2RAGErobustfunc(
        INV1final_array, INV2_array, noiselevel**2
    )

    # integer format
    if not integerformat == 0:
        MP2RAGERobustPhaseSensitive_array = np.round(
            4095 * (MP2RAGERobustPhaseSensitive_array + 0.5)
        )

    # save output
    fileOUT = os.path.split(fileUNI)
    output = nb.Nifti1Image(
        MP2RAGERobustPhaseSensitive_array, MP2RAGE_img.affine, MP2RAGE_img.header
    )
    nb.save(output, os.path.join(path_output, "n" + fileOUT[1]))
