# -*- coding: utf-8 -*-
"""Surface utilities."""

import datetime
import os
import subprocess
import sys

import numpy as np

from ..io.filename import get_filename

__all__ = ["mris_thickness", "mris_curvature", "inflate_surf_mesh"]


def mris_thickness(path, sub):
    """This function calculates the cortical thickness using freesurfer. Old files are
    moved to the freesurfer trash folder tagged with a date string as suffix.

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
    for _h in hemi:
        file_in = os.path.join(path, sub, "surf", f"{_h}.thickness")
        file_out = os.path.join(path_trash, f"{_h}.thickness_backup_{date}")
        os.rename(file_in, file_out)

        # get new thickness file
        command = "mris_thickness"
        command += f" -max {max_thickness}"
        command += f" {sub}"
        command += f" {_h}"
        command += f" {_h}.thickness"

        print("Execute: " + command)
        try:
            subprocess.run([command], shell=True, check=False)
        except subprocess.CalledProcessError:
            print("Execuation failed!")


def mris_curvature(file_in, path_output, a=10, dist=(10, 10), thresh=0.999):
    """This function calculates a curvature file for an input surface mesh using
    freesurfer. The input file needs to have a prefix which indicates the hemisphere of
    the surface mesh.

    Parameters
    ----------
    file_in : str
        Filename of input surface.
    path_output : str
        Path where output is written.
    a : int, optional
        Number of smoothing iterations. The default is 10.
    dist : tuple, optional
        Number of neighbours within distance in mm for curvature estimate, by default
        (10,10).
    thresh : float, optional
        Input threhold for curvature estimate, by default 0.999.

    Returns
    -------
    None.

    """
    # get hemi from filename
    hemi = os.path.splitext(os.path.basename(file_in))[0]
    if not hemi == "lh" and not hemi == "rh":
        sys.exit("Could not identify hemi from filename!")

    # calculate curvature file
    command = "mris_curvature"
    command += f" -a {a}"
    command += f" -distances {dist[0]} {dist[1]}"
    command += " -n -w"
    command += f" -thresh {thresh}"
    command += f" {file_in}"

    print("Execute: " + command)
    try:
        subprocess.run([command], shell=True, check=False)
    except subprocess.CalledProcessError:
        print("Execuation failed!")

    # rename mean curvature to curv
    os.rename(file_in + ".H", os.path.join(path_output, hemi + ".curv"))

    # delete gaussian curvature file and temporary input file
    os.remove(file_in + ".K")


def inflate_surf_mesh(file_in, file_out, n_iter):
    """The scripts takes a generated FreeSurfer surfaces and inflates it.

    Parameters
    ----------
    file_in : str
        Filename of input surface.
    file_out : str
        Filename of output surface.
    n_iter : int
        Number of inflating iterations.

    Returns
    -------
    None.

    """
    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # inflate surface
    try:
        subprocess.run(
            ["mris_inflate", "-n", str(n_iter), "-no-save-sulc", file_in, file_out],
            check=True,
        )
    except subprocess.CalledProcessError:
        sys.exit("Surface inflation failed!")


def match_vertex_number(vtx_white, vtx_pial, fac, ind_white, ind_pial):
    """This function takes arrays of white and pial surface vertices as input and
    matches the vertex numbers of both surfaces. The match is based on index lists which
    map the vertex indices to a common reference space. Removed vertices (not found in
    the corresponding other surface) are removed and the face array is updated. Index
    files are updated as well.

    Parameters
    ----------
    vtx_white : ndarray
        Vertex array of white surface.
    vtx_pial : ndarray
        Vertex array of pial surface.
    fac : ndarray
        Corresponding face array.
    ind_white : list
        Index list of white surface.
    ind_pial : list
        Index list of pial surface.

    Returns
    -------
    vtx_white : ndarray
        Updated vertex array of white surface.
    vtx_pial : ndarray
        Updated vertex array of pial surface.
    fac_new : ndarray
        Updated face array.
    ind_white : list
        updated index list.

    """
    # vertex indices in common reference space which are in neither of both
    # deformed surfaces
    ind_remove = list(set(ind_pial) - set(ind_white))
    ind_remove.extend(set(ind_white) - set(ind_pial))
    ind_remove = sorted(ind_remove, reverse=True)

    # sort vertices of white surface
    print("sort white surface vertices")

    i = 0
    c_white = np.zeros_like(ind_white)
    while i < len(ind_white):
        if np.any(ind_remove == ind_white[i]):
            c_white[i] = 1

        i += 1

    # sort vertices of pial surface
    print("sort pial surface vertices")

    i = 0
    c_pial = np.zeros_like(ind_pial)
    while i < len(ind_pial):
        if np.any(ind_remove == ind_pial[i]):
            c_pial[i] = 1

        i += 1

    # sort faces
    fac_old = fac.copy()
    fac_new = fac.copy()
    fac_outlier = np.zeros_like(fac)

    loop_status = 0
    loop_length = len(ind_remove)
    for i in range(loop_length):
        check_remove = np.where(ind_remove[i] == ind_white)[0]
        if len(check_remove):
            row, col = np.where(fac_old == check_remove)
            fac_outlier[row, col] = 1  # remember which faces to remove
            fac_temp = fac_new.copy()  # update face numbering
            fac_temp[fac_old > check_remove] = -1
            fac_temp[fac_temp != -1] = 0
            fac_new += fac_temp

        i += 1

        # print status
        counter = np.floor(i / loop_length * 100)
        if counter != loop_status:
            print("sort faces: " + str(counter) + " %")
            loop_status = counter

    # remove outliers in faces
    fac_outlier = np.sum(fac_outlier, 1)
    fac_new = fac_new[fac_outlier == 0]

    # remove outliers in vertices
    vtx_white = vtx_white[c_white == 0]
    vtx_pial = vtx_pial[c_pial == 0]

    # remove outliers in ind
    ind_white = ind_white[c_white == 0]

    return vtx_white, vtx_pial, fac_new, ind_white
