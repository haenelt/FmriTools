# -*- coding: utf-8 -*-

import datetime
import os


def get_thickness_fsurf(path, sub):
    """Get thickness fsurf.

    This function calculates the cortical thickness using freesurfer. Old files
    are moved to the freesurfer trash folder tagged with a date string as
    suffix.

    Parameters
    ----------
    path : str
        Path to the freesurfer segmentation folder.
    sub : str
        Name of the freesurfer segmentation folder.

    Returns
    -------
    None.

    """

    # parameters
    hemi = ["lh", "rh"]  # hemisphere prefix
    max_thickness = 10  # maximum thickness cut in mm

    # set subject environment here for mris_thickness
    os.environ["SUBJECTS_DIR"] = path

    # output folder (freesurfer trash folder)
    path_trash = os.path.join(path, sub, "trash")

    # get date string for moved files
    date = datetime.datetime.now().strftime("%Y%m%d%H%M")

    # move old surfaces and morphological files to trash folder
    for i in range(len(hemi)):
        file_in = os.path.join(path, sub, "surf", hemi[i] + ".thickness")
        file_out = os.path.join(path_trash, hemi[i] + ".thickness_backup_" + date)
        os.rename(file_in, file_out)

    # get new thickness file
    for i in range(len(hemi)):
        os.system(
            "mris_thickness"
            + " -max "
            + str(max_thickness)
            + " "
            + sub
            + " "
            + hemi[i]
            + " "
            + hemi[i]
            + ".thickness"
        )
