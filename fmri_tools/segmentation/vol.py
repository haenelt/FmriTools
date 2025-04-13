# -*- coding: utf-8 -*-
"""Volume utilities."""

import os

import nibabel as nb
import numpy as np
from scipy.ndimage import distance_transform_edt

from .. import execute_command
from ..io.filename import get_filename

__all__ = [
    "remove_bias_ants",
    "robust_combination",
    "include_pial_correction",
    "estimate_pv",
    "VDM",
]


def remove_bias_ants(file_in, file_out, save_bias=True, file_mask=None, strict=False):
    """Apply the ANTS N4 bias field correction.

    Parameters
    ----------
    file_in : str
        File name of input image.
    file_out : str
        File name of output image.
    save_bias : bool, optional
        Save bias field to disk.
    file_mask : str, optional
        File name of input mask.
    strict : bool, optional
        Strict mode.
    """
    # file name of saved bias field
    path, _, _ = get_filename(file_out)
    file_out_bias = os.path.join(path, "n4bias.nii")

    command = "N4BiasFieldCorrection"
    command += " -v 1 -d 3 -s 2 -r 0" if strict else " -d 3"
    command += f" --input-image {file_in}"
    if file_mask is not None:
        command += f" -x {file_mask}"
    if strict:
        command += " -c [50x50x50x50,0.0000001] -b [200]"
    if save_bias:
        command += f" --output [ {file_out}, {file_out_bias}]"
    else:
        command += f" --output {file_out}"

    # run
    execute_command(command)


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


def include_pial_correction(path, sub):
    """This function takes manual corrections from the file pial_edit.mgz which shall be
    located in the freesurfer mri folder and includes them in the brainmask. The manual
    corrected brainmask is saved as brain.finalsurfs.manedit.mgz. Note that the
    corrections in the pial_edit.mgz are done with brush value 256 and eraser value -1.

    Parameters
    ----------
    path : str
        Path to SUBJECTS_ID.
    sub : str
        Freesurfer subject name.
    """
    # open pial edit
    edit_img = nb.load(os.path.join(path, sub, "mri", "pial_edit.mgz"))
    edit_array = edit_img.get_fdata()

    # open brainmask
    brainmask_img = nb.load(os.path.join(path, sub, "mri", "brainmask.mgz"))
    brainmask_array = brainmask_img.get_fdata()

    # include manual edits to brainmask
    brainmask_array[edit_array == 256] = 255
    brainmask_array[edit_array == -1] = 1

    # save brainmask as brain.finalsurfs.manedit.mgz
    output = nb.Nifti1Image(brainmask_array, brainmask_img.affine, brainmask_img.header)
    nb.save(output, os.path.join(path, sub, "mri", "brain.finalsurfs.manedit.mgz"))


def estimate_pv(dir_out, subjects_dir, sub):
    """Use freesurfer mri_compute_volume_fractions to estimate partial volume
    contributions.

    Parameters
    ----------
    dir_out : str
        Output directory.
    subjects_dir : str
        Path to subject.
    sub : str
        Freesurfer subject name.
    """
    os.environ["SUBJECTS_DIR"] = subjects_dir

    command = "mri_compute_volume_fractions"
    command += f" --o {dir_out}/pv"
    command += " --regheader"
    command += f" {sub}"
    command += f" {subjects_dir}/{sub}/mri/brain.mgz"

    # run
    execute_command(command)


class VDM:
    """Vascular distance mapping (VDM) computes the euclidean distance of each voxel to
    its closest vessel. This approach has the advantage that a continuous metric can be
    used to analyze brain vasculature instead of a binary mask. The idea is taken from
    Bause et al. 2020, Mattern et al. (2020, 2021a, 2021b), Garcia-Garcia et al. (2023).

    Args:
        image: Nibabel image object containing the array of the binary vessel mask.

    References:
        - Bause, J., et al. Impact of prospective motion correction, distortion
          correction methods and large vein bias on the spatial accuracy of cortical
          laminar fMRI at 9.4 Tesla. NeuroImage (2020).
        - Mattern, H., et al. Vessel distance mapping. ESMRMB (2020).
        - Mattern, H., Vessel distance mapping of the aging subcortical venous
          vasculature. ESMRMB (2021a).
        - Mattern, H., et al. Vessel distance mapping for deep gray matter structures.
          ISMRM (2021b).
        - Garcia-Garcia, B., et al. Vessel distance mapping: A novel methodology for
          assessing vascular-induced cognitive resilience. NeuroImage (2023).

    """

    def __init__(self, image: nb.nifti1.Nifti1Image) -> None:
        self.arr = image.get_fdata()
        self.header = image.header
        self.affine = image.affine
        self.dim = image.header.get_zooms()
        self.distance_volume = None

    @property
    def transform(self):
        """Apply the distance transform."""
        self.distance_volume = distance_transform_edt(self.arr == 0, sampling=self.dim)
        return self.distance_volume

    def save(self, file_out: str) -> None:
        """Save VDM image to disk."""
        if self.distance_volume is None:
            print("Transform...")
            self.transform
        output = nb.Nifti1Image(self.distance_volume, self.affine, self.header)
        nb.save(output, file_out)

    @classmethod
    def from_file(cls, file_mask: str) -> "VDM":
        """Construct VDM object from file name."""
        image = nb.load(file_mask)
        return cls(image)
