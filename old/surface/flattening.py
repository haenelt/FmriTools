# -*- coding: utf-8 -*-

import datetime
import os
import shutil as sh

import numpy as np

from ..io.surf import label_as_patch

__all__ = ["label_flattening", "surface_flattening"]


def label_flattening(file_ref, file_label, path_output, cleanup=True):
    """Surface flattening.

    Uses the FreeSurfer mris_flatten function to flatten a dense patch of manually
    defined cortex. A label file is used to define the flattened region. Instead of
    smoothwm, we use white for surface flattening.

    Parameters
    ----------
    file_ref : str
        Reference surface file for flattening.
    file_label : str
        Label to be flattened saved as <hemi>.<namePATCH>.label.
    path_output : str
        Path where output is saved.
    cleanup : bool, optional
        Delete intermediate files. The default is True.

    Returns
    -------
    None.

    """

    # create temporary folder
    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = "".join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    path_temp = os.path.join(path_output, "tmp_" + tmp_string)

    # make temporary folder
    if not os.path.exists(path_temp):
        os.makedirs(path_temp)
    else:
        raise FileExistsError("Temporary folder already exists!")

    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # change to temporary folder
    cwd = os.getcwd()
    os.chdir(path_temp)

    # get hemisphere from file name
    hemi = os.path.basename(file_label)[:2]
    name_patch = "." + os.path.basename(file_label).split(".")[1]

    # convert label to patch file
    file_patch = os.path.join(path_temp, os.path.join(hemi + name_patch + ".patch.3d"))
    label_as_patch(file_ref, file_label, file_patch)

    # copy reference file and path into temporary folder
    sh.copy2(file_ref, os.path.join(path_temp, hemi + ".smoothwm"))

    # surface flattening
    w = 0  # write out the surface every number of iterations.
    s = 20  # size of neighbourhood to be used in the optimization
    n = 7  # number of vertices at each distance to be used in the optimization
    os.system(
        "mris_flatten"
        + " -w "
        + str(w)
        + " -distances "
        + str(s)
        + " "
        + str(n)
        + " "
        + file_patch
        + " "
        + hemi
        + name_patch
        + ".patch.flat"
    )

    # copy output
    sh.copy2(
        os.path.join(path_temp, hemi + name_patch + ".patch.flat"),
        os.path.join(path_output, hemi + name_patch + ".patch.flat"),
    )
    sh.copy2(
        os.path.join(path_temp, hemi + name_patch + ".patch.flat.out"),
        os.path.join(path_output, hemi + name_patch + ".patch.flat.out"),
    )

    # change to old folder
    os.chdir(cwd)

    # delete temporary files
    if cleanup:
        sh.rmtree(path_temp, ignore_errors=True)


def surface_flattening(file_ref, file_patch, path_output, cleanup=True):
    """Surface flattening.

    Uses the FreeSurfer mris_flatten function to flatten a dense patch of
    manually defined cortex. The manually cutted patch should have the following
    file name <hemi>.<namePATCH>.patch.3d. Instead of smoothwm, we use white for
    surface flattening.

    Parameters
    ----------
    file_ref : str
        Reference surface file for flattening.
    file_patch : str
        Patch to be flattened saved as <hemi>.<namePATCH>.patch.3d.
    path_output : str
        Path where output is saved.
    cleanup : bool, optional
        Delete intermediate files. The default is True.

    Returns
    -------
    None.

    """

    # create temporary folder
    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = "".join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    path_temp = os.path.join(path_output, "tmp_" + tmp_string)

    # make temporary folder
    if not os.path.exists(path_temp):
        os.makedirs(path_temp)
    else:
        raise FileExistsError("Temporary folder already exists!")

    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # change to temporary folder
    cwd = os.getcwd()
    os.chdir(path_temp)

    # divide patch basename
    hemi = os.path.splitext(
        os.path.splitext(os.path.splitext(os.path.basename(file_patch))[0])[0]
    )[0]
    name_patch = os.path.splitext(
        os.path.splitext(os.path.splitext(os.path.basename(file_patch))[0])[0]
    )[1]

    # copy reference file and path into temporary folder
    sh.copy2(file_ref, os.path.join(path_temp, hemi + ".smoothwm"))
    sh.copy2(file_patch, os.path.basename(file_patch))

    # surface flattening
    w = 0  # write out the surface every number of iterations.
    s = 20  # size of neighbourhood to be used in the optimization
    n = 7  # number of vertices at each distance to be used in the optimization
    os.system(
        "mris_flatten"
        + " -w "
        + str(w)
        + " -distances "
        + str(s)
        + " "
        + str(n)
        + " "
        + hemi
        + name_patch
        + ".patch.3d"
        + " "
        + hemi
        + name_patch
        + ".patch.flat"
    )

    # copy output
    sh.copy2(
        os.path.join(path_temp, hemi + name_patch + ".patch.flat"),
        os.path.join(path_output, hemi + name_patch + ".patch.flat"),
    )
    sh.copy2(
        os.path.join(path_temp, hemi + name_patch + ".patch.flat.out"),
        os.path.join(path_output, hemi + name_patch + ".patch.flat.out"),
    )

    # change to old folder
    os.chdir(cwd)

    # delete temporary files
    if cleanup:
        sh.rmtree(path_temp, ignore_errors=True)
