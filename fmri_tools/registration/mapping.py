# -*- coding: utf-8 -*-
"""Map to surface."""

import datetime
import os
import shutil as sh
import subprocess
import sys

import nibabel as nb
import numpy as np
from nibabel.freesurfer.io import read_geometry, read_morph_data, write_morph_data
from scipy.interpolate import griddata
from scipy.ndimage.filters import gaussian_filter

from ..io.affine import vox2ras_tkr
from ..io.filename import get_filename
from ..io.surf import read_mgh, write_mgh
from ..io.vol import mri_convert
from ..registration.transform import apply_affine_chunked
from ..utils.interpolation import linear_interpolation3d, nn_interpolation3d

__all__ = [
    "mri_vol2surf",
    "map2surface",
    "map_timeseries",
    "map2grid",
    "map2stack",
    "morph2dense",
]

# linear and nearest neighbor sampling functions
_sampler = {"linear": linear_interpolation3d, "nearest": nn_interpolation3d}


def mri_vol2surf(file_in, file_out, sub, interp_method="nearest"):
    """Use freesurfer mri_vol2surf to sample volume data onto a surface mesh.

    Parameters
    ----------
    file_in : str
        File name of input volume file.
    file_out : str
        File name of output overlay file.
    sub : str
        Name of freesurfer subject.
    interp_method : str, optional
        Interpolation method (nearest or trilinear), by default "nearest"
    """
    # get hemisphere
    _, _, hemi = get_filename(file_in)
    hemi = hemi.replace(".", "")
    if hemi not in ["lh", "rh"]:
        raise ValueError("No hemisphere specified in file name!")

    command = "mri_vol2surf"
    command += f" --hemi {hemi}"
    command += f" --interp {interp_method}"
    command += f" --o {file_out}"
    command += " --out_type mgh"
    command += f" --regheader {sub}"
    command += " --projdist 0.000"
    command += f" --mov {file_in}"
    command += " --surf source"

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")


def map2surface(
    input_surf,
    input_vol,
    write_output=False,
    path_output="",
    interp_method="nearest",
    input_surf_target=None,
    input_ind=None,
    cleanup=True,
):
    """This function samples data from the input volume to the input surface and
    optionally maps those values to a target surface if an index file is given.

    Parameters
    ----------
    input_surf : str
        Surface mesh onto which volume data is sampled.
    input_vol : str
        Volume from which data is sampled.
    write_output : bool, optional
        Write sampled data as MGH file. The default is False.
    path_output : str, optional
        Path where to save output. The default is "".
    interp_method : str, optional
        Interpolation method (nearest or trilinear). The default is "nearest".
    input_surf_target : str, optional
        Target surface (only necessary if index file is given). The default is
        None.
    input_ind : str, optional
        Textfile with mapping of vertex indices to target space. The default is
        None.
    cleanup : bool, optional
        Remove intermediate files. The default is True.

    Raises
    ------
    FileExistsError
        If the temporary folder already exists.

    Returns
    -------
    arr_sampled : ndarray
        Image array.
    affine_sampled : ndarray
        Affine transformation matrix.
    header_sampled : MGHHeader
        Image header.

    """
    # clean everything if no output is written
    if not write_output:
        path_output, _, _ = get_filename(input_vol)
        # cleanup = True

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
        os.makedirs(path_output)

    # mimic freesurfer folder structure (with some additional folder for
    # intermediate files)
    path_sub = os.path.join(path_output, sub)
    path_mri = os.path.join(path_sub, "mri")
    path_surf = os.path.join(path_sub, "surf")

    if not os.path.exists(path_sub):
        os.makedirs(path_sub)
    else:
        raise FileExistsError("Temporary folder already exists!")

    os.makedirs(path_mri)
    os.makedirs(path_surf)

    # get filenames
    _, name_vol, ext_vol = get_filename(input_vol)
    _, hemi, name_surf = get_filename(input_surf)
    name_surf = name_surf.replace(".", "")

    # check filename
    if not hemi == "lh" and not hemi == "rh":
        sys.exit("Could not identify hemi from filename!")

    # copy input volume as orig.mgz to mimic freesurfer folder
    if ext_vol != ".mgz":
        mri_convert(input_vol, os.path.join(path_mri, name_vol + ".mgz"))
        os.rename(
            os.path.join(path_mri, name_vol + ".mgz"),
            os.path.join(path_mri, "orig.mgz"),
        )
    else:
        sh.copyfile(input_vol, os.path.join(path_mri, "orig.mgz"))

    # copy input surface to mimic freesurfer folder
    sh.copyfile(input_surf, os.path.join(path_surf, hemi + ".source"))

    # filename of sampled data
    file_sampled = os.path.join(path_surf, hemi + "." + "sampled.mgh")

    # mri_vol2surf
    mri_vol2surf(input_vol, file_sampled, sub, interp_method=interp_method)

    # load data
    arr_sampled, affine_sampled, header_sampled = read_mgh(file_sampled)

    # map on separate mesh
    if input_ind:
        # load data
        ind_target = np.loadtxt(input_ind, dtype=int)
        vtx_target, _ = read_geometry(input_surf_target)

        # read sampled morph data
        arr_tmp = arr_sampled.copy()

        # update header
        header_sampled["dims"][0] = len(vtx_target)
        header_sampled["Mdc"] = np.eye(3)

        # sample array in target space
        arr_sampled = np.zeros(len(vtx_target))
        arr_sampled[ind_target] = arr_tmp

    if write_output:
        file_out = os.path.join(path_output, hemi + "." + name_vol + "_" + name_surf)

        if input_ind:
            file_out += "_trans.mgh"
        else:
            file_out += ".mgh"

        write_mgh(file_out, arr_sampled, affine_sampled, header_sampled)

        # delete intermediate files
    if cleanup:
        sh.rmtree(path_sub, ignore_errors=True)
        if not len(os.listdir(path_output)):
            sh.rmtree(path_output)

    return arr_sampled, affine_sampled, header_sampled


def map_timeseries(vtx, arr_timeseries, dims, ds, interpolation="linear"):
    """Map time series data onto a surface mesh. A 2D array is returned which contains
    vertex-wise sampled data for each time point in separate columns. All vertices
    outside the time series volume are set to nan.

    Parameters
    ----------
    vtx : np.ndarray, shape=(N,3)
        Vertex-wise array.
    arr_timeseries : np.ndarray, shape=(X,Y,Z,T)
        4D array of fMRI time series.
    dims : tuple
        Tuple containing volume dimensions in x-, y- and z-direction.
    ds : tuple
        Tuple containing voxel sizes in x-, y- and z-direction.
    interpolation : str, optional (linear | nearest)
        Interpolation method (linear or nearest neighbor interpolation).

    Returns
    -------
    arr_sampled : np.ndarray, shape=(N,T)
        Vertex-wise sampled time series.

    """
    nx, ny, nz, nt = np.shape(arr_timeseries)
    _, ras2vox = vox2ras_tkr(dims, ds)
    vtx_vox = apply_affine_chunked(ras2vox, vtx)

    # exclude nans and vertices outside of the volume
    mask = np.ones(len(vtx), dtype=bool)
    for i, n in enumerate((nx, ny, nz)):
        mask[np.isnan(vtx_vox[:, i])] = 0
        mask[vtx_vox[:, i] < 0] = 0
        mask[vtx_vox[:, i] > n - 1] = 0
    vtx_vox = vtx_vox[mask == 1, :]

    # initialize resulting array
    arr_sampled = np.empty((len(vtx), nt))
    arr_sampled[:] = np.nan

    for i in range(nt):
        arr_sampled[mask == 1, i] = _sampler[interpolation](
            vtx_vox[:, 0], vtx_vox[:, 1], vtx_vox[:, 2], arr_timeseries[:, :, :, i]
        )

    return arr_sampled


def map2grid(
    file_grid,
    file_input,
    sigma,
    path_output="",
    basename_output="",
    binary=False,
    overwrite=True,
):
    """This script allows you to sample indexed morphological data onto the regular
    grid. Optional, a gaussian filter can be applied to the output image.

    Parameters
    ----------
    file_grid : str
        Filename of grid coordinate mapping.
    file_input : str
        Filename of morphological data or *.mgh data.
    sigma : float
        Standard deviation of Gaussian kernel.
    path_output : str, optional
        Path where output is saved. The default is "".
    basename_output : str, optional
        Basename of written output file. The default is "".
    binary : bool, optional
        Threshold output grid (for curvature file). The default is False.
    overwrite : bool, optional
        Write output to file. The default is True.

    Returns
    -------
    grid_array : ndarray
        File mapped onto array.

    """
    # load data
    grid_img = nb.load(file_grid)
    grid_array = grid_img.get_fdata()
    if os.path.splitext(file_input)[1] == ".mgh":
        morph = nb.load(file_input).get_fdata()
    else:
        morph = read_morph_data(file_input)

    # sample data onto grid
    for i in range(np.size(grid_array, 0)):
        for j in range(np.size(grid_array, 1)):
            if grid_array[i, j] != 0:
                grid_array[i, j] = morph[grid_array[i, j].astype(int)]

    # gaussian filter (opt)
    if sigma != 0:
        order = 0
        mode = "reflect"
        truncate = 4.0
        grid_array = gaussian_filter(
            grid_array, sigma=sigma, order=order, mode=mode, truncate=truncate
        )

    # binary mode (opt)
    if binary is True:
        grid_array[grid_array > 0] = 1
        grid_array[grid_array != 1] = -1

    # write output data
    if overwrite:
        # make output folder
        if not os.path.exists(path_output):
            os.mkdir(path_output)

        if sigma == 0 and binary is True:
            filenameOUT = os.path.join(
                path_output, basename_output + "_grid_binary.nii"
            )
        elif sigma == 0 and binary is False:
            filenameOUT = os.path.join(path_output, basename_output + "_grid.nii")
        elif sigma != 0 and binary is True:
            filenameOUT = os.path.join(
                path_output,
                basename_output + "_sigma" + str(sigma) + "_grid_binary.nii",
            )
        else:
            filenameOUT = os.path.join(
                path_output, basename_output + "_sigma" + str(sigma) + "_grid.nii"
            )

        output = nb.Nifti1Image(grid_array, grid_img.affine, grid_img.header)
        nb.save(output, filenameOUT)

    return grid_array


def map2stack(file_data, file_grid, sigma, path_output):
    """This function allows you to sample surface data to a patch defined on a regular
    grid. If multiple data files are given in a list, all grids are stacked together.

    Parameters
    ----------
    file_data : str
        Filename list of data.
    file_grid : str
        Filename of grid coordinate mapping.
    sigma : float
        Standard deviation of Gaussian kernel.
    path_output : str
        Path where output is saved.

    Returns
    -------
    None.

    """
    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # load data
    grid_img = nb.load(file_grid)
    grid_array = grid_img.get_fdata()

    # dim
    x = grid_img.header["dim"][1]
    y = grid_img.header["dim"][2]
    z = len(file_data)

    # sample data onto grid
    stack_array = np.zeros((x, y, z))
    for i in range(z):
        # load data
        data_img = nb.load(file_data[i])
        data_array = data_img.get_fdata()

        # map to stack
        for j in range(x):
            for k in range(y):
                if grid_array[j, k] != 0:
                    stack_array[j, k, i] = data_array[grid_array[j, k].astype(int)]

        # gaussian filter (opt)
        if sigma != 0:
            order = 0
            mode = "reflect"
            truncate = 4.0
            for j in range(np.size(stack_array, 2)):
                stack_array[:, :, j] = gaussian_filter(
                    stack_array[:, :, j],
                    sigma=sigma,
                    order=order,
                    mode=mode,
                    truncate=truncate,
                )

    # write output data
    filename_out = os.path.join(
        path_output,
        os.path.splitext(os.path.basename(file_data[0]))[0]
        + "_sigma"
        + str(sigma)
        + "_grid.nii",
    )
    output = nb.Nifti1Image(stack_array, grid_img.affine, grid_img.header)
    nb.save(output, filename_out)


def morph2dense(source_sphere, target_sphere, input_morph, path_output):
    """This function maps a morphological file from a source to a target surface.

    Parameters
    ----------
    source_sphere : str
        Source surface.
    target_sphere : str
        Target surface.
    input_morph : str
        Morphological input file.
    path_output : str
        Path where output is saved.

    Returns
    -------
    None.

    """
    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # transform morphological data to dense surfaces
    pts_sphere_dense, _ = read_geometry(target_sphere)
    pts_sphere, _ = read_geometry(source_sphere)

    # get morphological data
    morph = read_morph_data(input_morph)

    # do the transformation
    method = "nearest"
    morph_dense = griddata(pts_sphere, morph, pts_sphere_dense, method)

    # write dense morphological data
    write_morph_data(
        os.path.join(path_output, os.path.basename(input_morph)), morph_dense
    )
