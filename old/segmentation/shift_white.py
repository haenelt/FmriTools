# -*- coding: utf-8 -*-

import datetime
import os
import shutil

from nipype.interfaces.freesurfer import MRIsExpand, MRIsInflate

from ..segmentation.surf import mris_curvature


def shift_white(path, sub, w_shift=-0.5):
    """Shift white.

    This function shifts the white matter surface inwards (if negative values
    are applied). This step is included since the white matter surface is biased
    and lies somewhat within grey matter for MP2RAGE segmentations. A shift of
    -0.5 mm seems to be a reasonable value. After applying the shift, some
    measurements are updated including hemi.inflated, hemi.sulc and hemi.curv.
    Old files are moved to the freesurfer trash folder tagged with a date string
    as suffix. The curvature file is computed as mean curvature (not Gaussian
    curvature).

    Parameters
    ----------
    path : str
        Path to the freesurfer segmentation folder.
    sub : str
        Name of the freesurfer segmentation folder.
    w_shift : float, optional
        Inward shift (negative) of white surface in mm. The default is -0.5.

    Returns
    -------
    None.

    """

    # parameters
    hemi = ["lh", "rh"]  # hemisphere prefix
    n_inflate = 50  # number of iterations for surface inflation
    curv_distance = (
        10,
        10,
    )  # number of neighbours within distance in mm for curvature estimate
    curv_threshold = 0.999  # input threhold for curvature estimate
    curv_average = 5  # number of smoothing iterations for curvature estimate

    # output folder (freesurfer trash folder)
    path_trash = os.path.join(path, sub, "trash")

    # filenames to be moved
    filename = ["white", "inflated", "curv", "sulc"]

    # get date string for moved files
    date = datetime.datetime.now().strftime("%Y%m%d%H%M")

    # move old surfaces and morphological files to trash folder
    for i in range(len(hemi)):
        for j in range(len(filename)):
            file_in = os.path.join(path, sub, "surf", hemi[i] + "." + filename[j])
            file_out = os.path.join(
                path_trash, hemi[i] + "." + filename[j] + "_backup_" + date
            )
            os.rename(file_in, file_out)

    # shift white matter surface inwards
    for i in range(len(hemi)):
        mris_expand = MRIsExpand()
        mris_expand.inputs.in_file = os.path.join(
            path_trash, hemi[i] + ".white_backup_" + date
        )
        mris_expand.inputs.thickness = False
        mris_expand.inputs.distance = w_shift
        mris_expand.inputs.out_name = os.path.join(
            path, sub, "surf", hemi[i] + ".white"
        )
        mris_expand.run()

    # get new inflation and sulc file
    for i in range(len(hemi)):
        inflate = MRIsInflate()
        inflate.inputs.in_file = os.path.join(path, sub, "surf", hemi[i] + ".white")
        inflate.inputs.args = "-n " + str(n_inflate)
        inflate.inputs.out_file = os.path.join(path, sub, "surf", hemi[i] + ".inflated")
        inflate.inputs.out_sulc = os.path.join(path, sub, "surf", hemi[i] + ".sulc")
        inflate.run()

    # get new curvature file
    for _h in hemi:
        # copy the input file to not overwrite existing files
        file_in = os.path.join(path, sub, "surf", _h + ".white")
        file_out = os.path.join(path, sub, "surf", _h + ".temp")
        shutil.copy(file_in, file_out)

        mris_curvature(
            os.path.join(path, sub, "surf", _h + ".temp"),
            os.path.join(path, sub, "surf"),
            curv_average,
            curv_distance,
            curv_threshold,
        )

        # delete intermediate files afterwards
        os.remove(os.path.join(path, sub, "surf", _h + ".temp"))
