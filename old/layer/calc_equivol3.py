# -*- coding: utf-8 -*-

import os

import nibabel as nb
from nighres.laminar import volumetric_layering
from nighres.surface import probability_to_levelset

from ..io.vol import mri_convert


def calc_equivol3(path, sub, n_layers, path_output):
    """Calc equivol 3.

    This script computes the volumetric layering from the FreeSurfer ribbon
    mask.

    Parameters
    ----------
    path : str
        Path to the freesurfer segmentation folder.
    sub : str
        Name of the freesurfer segmentation folder.
    n_layers : int
        Number of layers.
    path_output : str
        Path where output is saved.

    Returns
    -------
    None.

    """

    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    if not os.path.exists(os.path.join(path_output, "binary")):
        os.mkdir(os.path.join(path_output, "binary"))

    if not os.path.exists(os.path.join(path_output, "levelset")):
        os.mkdir(os.path.join(path_output, "levelset"))

    # convert ribbon to nifti
    mri_convert(
        os.path.join(path, sub, "mri", "ribbon.mgz"),
        os.path.join(path, sub, "mri", "ribbon.nii.gz"),
    )

    # get GM/WM and GM/CSF boundaries
    ribbon_img = nb.load(os.path.join(path, sub, "mri", "ribbon.nii.gz"))
    ribbon_array = ribbon_img.get_fdata()

    # layering using nighres
    hemi = ["lh", "rh"]
    wm = [2, 41]
    gm = [3, 42]
    for i in range(len(hemi)):
        wm_array = ribbon_array.copy()
        wm_array[wm_array == wm[i]] = 1
        wm_array[wm_array != 1] = 0

        csf_array = ribbon_array.copy()
        csf_array[csf_array == gm[i]] = 1
        csf_array[csf_array == wm[i]] = 1
        csf_array[csf_array != 1] = 0

        output = nb.Nifti1Image(wm_array, ribbon_img.affine, ribbon_img.header)
        nb.save(output, os.path.join(path_output, "binary", hemi[i] + "_wm.nii"))

        output = nb.Nifti1Image(csf_array, ribbon_img.affine, ribbon_img.header)
        nb.save(output, os.path.join(path_output, "binary", hemi[i] + "_csf.nii"))

        # probability to levelset
        probability_to_levelset(
            os.path.join(path_output, "binary", hemi[i] + "_wm.nii"),
            save_data=True,
            output_dir=os.path.join(path_output, "levelset"),
            file_name=hemi[i] + "_wm",
        )

        # probability to levelset
        probability_to_levelset(
            os.path.join(path_output, "binary", hemi[i] + "_csf.nii"),
            save_data=True,
            output_dir=os.path.join(path_output, "levelset"),
            file_name=hemi[i] + "_csf",
        )

        # layering
        suffix = "_p2l-surf"
        volumetric_layering(
            os.path.join(path_output, "levelset", hemi[i] + "_wm" + suffix + ".nii.gz"),
            os.path.join(
                path_output, "levelset", hemi[i] + "_csf" + suffix + ".nii.gz"
            ),
            n_layers=n_layers,
            topology_lut_dir=None,
            save_data=True,
            output_dir=path_output,
            file_name=hemi[i],
        )
