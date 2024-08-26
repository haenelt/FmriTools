# -*- coding: utf-8 -*-
"""Map to surface."""

import datetime
import os
import shutil as sh
import subprocess
import sys

import nibabel as nb
import numpy as np
from nibabel.freesurfer.io import (read_geometry, read_morph_data,
                                   write_geometry, write_morph_data)
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from sh import gunzip

from ..io.affine import apply_affine_chunked, read_vox2ras_tkr, vox2ras_tkr
from ..io.filename import get_filename
from ..io.surf import read_mgh, write_mgh
from ..io.vol import mri_convert
from ..preprocessing.fieldmap import fugue_fsl, prepare_fieldmap_fsl
from ..registration.cmap import generate_coordinate_mapping
from ..registration.fsl import flirt
from ..segmentation.skullstrip import skullstrip_epi
from ..surface.mesh import Mesh
from ..surface.smooth import mris_smooth
from ..utils.interpolation import linear_interpolation3d, nn_interpolation3d
from ..utils.roi import erode_fsl
from .fsl import apply_flirt, apply_fugue

__all__ = [
    "mri_vol2surf",
    "map2surface",
    "map_timeseries",
    "map2grid",
    "map2stack",
    "morph2dense",
    "mesh_sampling",
    "deform_surface",
    "remove_vertex_outliers",
    "apply_fieldmap",
]

# linear and nearest neighbor sampling functions
_sampler = {"linear": linear_interpolation3d, "nearest": nn_interpolation3d}


def mri_vol2surf(
    file_in, file_out, subjects_dir, sub, hemi, surf="white", interp_method="nearest"
):
    """Use freesurfer mri_vol2surf to sample volume data onto a surface mesh.

    Parameters
    ----------
    file_in : str
        File name of input volume file.
    file_out : str
        File name of output overlay file.
    subjects_dir : str
        Path to subject.
    sub : str
        Name of freesurfer subject.
    hemi : str
        Hemisphere.
    surf : str
        Surface name.
    interp_method : str, optional
        Interpolation method (nearest or trilinear), by default "nearest"
    """
    os.environ["SUBJECTS_DIR"] = subjects_dir

    command = "mri_vol2surf"
    command += f" --hemi {hemi}"
    command += f" --interp {interp_method}"
    command += f" --o {file_out}"
    command += " --out_type mgh"
    command += f" --regheader {sub}"
    command += " --projdist 0.000"
    command += f" --mov {file_in}"
    command += f" --surf {surf}"

    print("Execute: " + command)
    try:
        subprocess.run([command], shell=True, check=False)
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


def mesh_sampling(
    surf_in,
    vol_in,
    write_output=False,
    path_output="",
    source2target_in="",
    interp_method="nearest",
    cleanup=True,
):
    """This function samples data onto a surface mesh. Optionally, a coordinate mapping
    can be applied to transform the surface mesh to the space of the input volume.

    Parameters
    ----------
    surf_in : str
        Filename of input surface mesh.
    vol_in : str
        Filename of input volume from which data is sampled.
    write_output : bool, optional
        Write sampled data as MGH file. The default is False.
    path_output : str, optional
        Path where output is written. The default is "".
    source2target_in : str, optional
        Source to target coordinate mapping. The default is "".
    interp_method : str, optional
        Interpolation method for surface sampling. Possible arguments are
        nearest and trilinear. The default is "nearest".
    cleanup : bool, optional
        Remove intermediate files. The default is True.

    Returns
    -------
    arr : ndarray
        Image array.
    affine : ndarray
        Affine transformation matrix.
    header : MGHHeader
        Image header.

    """
    # clean everything if no output is written
    if not write_output:
        path_output, _, _ = get_filename(vol_in)
        # cleanup = True

    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = "".join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    path_tmp = os.path.join(path_output, "tmp_" + tmp_string)
    if not os.path.exists(path_tmp):
        os.makedirs(path_tmp)
    else:
        raise FileExistsError("Temporary folder already exists!")

    # get filenames
    _, hemi, name_mesh = get_filename(surf_in)
    name_mesh = name_mesh.replace(".", "")
    _, name_vol, ext_vol = get_filename(vol_in)

    # check filename
    if not hemi == "lh" and not hemi == "rh":
        sys.exit("Could not identify hemi from filename!")

    # copy temporary vol
    file_vol = os.path.join(path_tmp, name_vol + ext_vol)
    sh.copy(vol_in, file_vol)

    # unzip if necessary
    if file_vol[-3:] == ".gz":
        gunzip(file_vol)
        file_vol = file_vol[:-3]

    # get cmap
    if source2target_in:
        _, name_s2t, ext_s2t = get_filename(source2target_in)
        file_s2t = os.path.join(path_tmp, name_s2t + ext_s2t)
        sh.copy(source2target_in, file_s2t)
    else:
        name_s2t = "cmap_t2t"
        ext_s2t = ".nii"
        generate_coordinate_mapping(
            vol_in,
            pad=0,
            path_output=path_tmp,
            suffix="t2t",
            time=False,
            write_output=True,
        )
        file_s2t = os.path.join(path_tmp, name_s2t + ext_s2t)

    if file_s2t[-3:] == ".gz":
        gunzip(file_s2t)
        file_s2t = file_s2t[:-3]

    # deform mesh
    deform_surface(
        input_surf=surf_in,
        input_orig=file_s2t,
        input_deform=file_s2t,
        input_target=file_vol,
        path_output=path_tmp,
        input_mask=None,
        interp_method="trilinear",
        smooth_iter=0,
        flip_faces=False,
        cleanup=True,
    )

    # do mapping
    file_def = os.path.join(path_tmp, hemi + "." + name_mesh + "_def")
    arr, affine, header = map2surface(
        input_surf=file_def,
        input_vol=file_vol,
        write_output=False,
        path_output=path_tmp,
        interp_method=interp_method,
        input_surf_target=None,
        input_ind=None,
        cleanup=True,
    )

    if write_output:
        _, name_vol, _ = get_filename(file_vol)
        file_out = os.path.join(
            path_output, hemi + "." + name_vol + "_" + name_mesh + ".mgh"
        )
        write_mgh(file_out, arr, affine, header)

    # delete intermediate files
    if cleanup:
        sh.rmtree(path_tmp, ignore_errors=True)
        if not len(os.listdir(path_output)):
            sh.rmtree(path_output)

    return arr, affine, header


def deform_surface(
    input_surf,
    input_orig,
    input_deform,
    input_target,
    path_output,
    input_mask=None,
    interp_method="nearest",
    smooth_iter=0,
    flip_faces=False,
    cleanup=True,
):
    """This function deforms a surface mesh in freesurfer convention using a coordinate
    map containing voxel coordinates. The computation takes quite a while because in the
    case of removed vertices, i.e. if a mask is given as input, the remaining faces are
    reindexed.

    Parameters
    ----------
    input_surf : str
        Surface mesh to be transformed.
    input_orig : str
        Freesurfer orig.mgz.
    input_deform : str
        Deformation (coordinate mapping).
    input_target : str
        Target volume.
    path_output : str
        Path where to save output.
    input_mask : str, optional
        Mask volume. The default is None.
    interp_method : str, optional
        Interpolation method (nearest or trilinear). The default is "nearest".
    smooth_iter : int, optional
        Number of smoothing iterations applied to final image (if set > 0). The
        default is 0.
    flip_faces : bool, optional
        Reverse normal direction of mesh. The default is False.
    cleanup : bool, optional
        Remove intermediate files. The default is True.
    """
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
    _, _, ext_orig = get_filename(input_orig)
    _, hemi, name_surf = get_filename(input_surf)
    name_surf = name_surf.replace(".", "")

    # check filename
    if not hemi == "lh" and not hemi == "rh":
        sys.exit("Could not identify hemi from filename!")

    # copy orig, cmap and input surface to mimicked freesurfer folders
    sh.copyfile(input_surf, os.path.join(path_surf, hemi + ".source"))
    if ext_orig != ".mgz":
        mri_convert(input_orig, os.path.join(path_mri, "orig.mgz"))
    else:
        sh.copyfile(input_orig, os.path.join(path_mri, "orig.mgz"))

    # read surface geometry
    vtx, fac = read_geometry(input_surf)

    # get affine vox2ras-tkr transformation to target volume
    vox2ras_tkr, _ = read_vox2ras_tkr(input_target)

    # divide coordinate mapping into its x, y and z components
    cmap_img = nb.load(input_deform)
    cmap_img.header["dim"][0] = 3
    cmap_img.header["dim"][4] = 1

    # apply vox2ras transformation to coordinate mappings
    cmap_array = cmap_img.get_fdata()
    cmap_array = apply_affine_chunked(vox2ras_tkr, cmap_array)

    components = ["x", "y", "z"]
    vtx_new = np.zeros([len(vtx), 3])
    for i, _ in enumerate(components):
        file_temp = os.path.join(path_mri, components[i] + "_deform.nii")
        file_sampled = os.path.join(
            path_surf, hemi + "." + components[i] + "_sampled.mgh"
        )

        # get target volume
        temp_array = cmap_array[:, :, :, i]
        temp_img = nb.Nifti1Image(temp_array, cmap_img.affine, cmap_img.header)
        nb.save(temp_img, file_temp)

        # mri_vol2surf
        mri_vol2surf(file_temp, file_sampled, sub, interp_method=interp_method)

        data_img = nb.load(file_sampled)
        vtx_new[:, i] = np.squeeze(data_img.get_fdata())

    if input_mask:
        file_background = os.path.join(path_surf, hemi + ".background.mgh")

        # mri_vol2surf (background)
        mri_vol2surf(input_mask, file_background, sub, interp_method="nearest")

        # get new indices
        background_list = nb.load(file_background).get_fdata()
        background_list = np.squeeze(background_list).astype(int)

        # only keep vertex indices within the slab
        ind_keep = np.arange(len(vtx))
        ind_keep = ind_keep[background_list != 0]

        vtx_new, fac_new, ind_keep = Mesh(vtx_new, fac).remove_vertices(
            ind_keep, create_ind=True
        )

        # save index mapping between original and transformed surface
        np.savetxt(
            os.path.join(path_output, hemi + "." + name_surf + "_ind.txt"),
            ind_keep,
            fmt="%d",
        )
    else:
        fac_new = fac

    # flip faces
    if flip_faces:
        fac_new = np.flip(fac_new, axis=1)

    # write new surface
    file_out = os.path.join(path_output, hemi + "." + name_surf + "_def")
    write_geometry(file_out, vtx_new, fac_new)

    # smooth surface
    if smooth_iter:
        mris_smooth(file_out, file_out + "_smooth", smooth_iter)

    # delete intermediate files
    if cleanup:
        sh.rmtree(path_sub, ignore_errors=True)


def remove_vertex_outliers(input_surf, input_ind, n=5, overwrite=True):
    """This function removes outlier vertices from a deformed surface mesh. Due to
    interpolation of the deformation fields, edge effects can move some vertices to the
    edge of the image volume. These are removed by comparing each vertex to the
    geometric center of the whole mesh and setting a global threshold. The threshold is
    defined as multiple (n) of the standard deviation of the distance distribution. This
    is a very crude method to delete outlier vertices.

    Parameters
    ----------
    input_surf : str
        Path to the input surface mesh.
    input_ind : str
        Path to the corresponding index list (with .txt extension).
    n : TYPE, float
        Threshold parameter. The default is 5.
    overwrite : bool, optional
        Overwrite input surface. The default is True.
    """
    # load geometry
    vtx, fac = read_geometry(input_surf)

    # load index file
    ind = np.loadtxt(input_ind)

    # get geometry center
    x_mean = np.sum(vtx[:, 0]) / len(vtx[:, 0])
    y_mean = np.sum(vtx[:, 1]) / len(vtx[:, 1])
    z_mean = np.sum(vtx[:, 2]) / len(vtx[:, 2])

    # euclidean distance to geometric center
    vtx_dist = np.zeros_like(vtx)
    vtx_dist[:, 0] = (vtx[:, 0] - x_mean) ** 2
    vtx_dist[:, 1] = (vtx[:, 1] - y_mean) ** 2
    vtx_dist[:, 2] = (vtx[:, 2] - z_mean) ** 2
    vtx_dist = np.sqrt(np.sum(vtx_dist, 1))

    # mean and std distance
    vtx_dist_mean = np.mean(vtx_dist)
    vtx_dist_std = np.std(vtx_dist)

    # distance threshold
    vtx_dist_threshold = vtx_dist_mean + n * vtx_dist_std

    # sort faces
    fac_old = fac.copy()
    fac_outlier = np.zeros_like(fac)
    n_outlier = np.zeros(len(vtx))
    c_step = 0
    n_step = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    for i in range(len(vtx)):
        if vtx_dist[i] > vtx_dist_threshold:
            row, col = np.where(fac_old == i)
            fac_outlier[row, col] = 1  # remember which faces to remove
            n_outlier[i] = 1  # remember which vertices to remove
            fac_temp = fac.copy()  # update face numbering
            fac_temp[fac_old >= i] = -1
            fac_temp[fac_temp != -1] = 0
            fac += fac_temp

        # print status
        counter = np.floor(i / len(vtx) * 100).astype(int)
        if counter == n_step[c_step]:
            print("remove outliers: " + str(counter) + " %")
            c_step += 1

    # remove outlier faces
    fac_outlier = np.sum(fac_outlier, 1)
    fac = fac[fac_outlier == 0]

    # remove outliers in vertex and ind
    vtx = vtx[n_outlier == 0]
    ind = ind[n_outlier == 0]

    # write output
    if overwrite:
        write_geometry(input_surf, vtx, fac)
        np.savetxt(input_ind, ind, fmt="%d")
    else:
        path_output = os.path.dirname(input_surf)
        name_output = os.path.basename(input_surf)
        write_geometry(os.path.join(path_output, name_output + "_out"), vtx, fac)

        path_output = os.path.dirname(input_ind)
        name_output = os.path.splitext(os.path.basename(input_ind))[0]
        np.savetxt(os.path.join(path_output, name_output + "_out.txt"), ind, fmt="%d")


def apply_fieldmap(
    file_fmap_magn,
    file_fmap_phase,
    file_epi,
    file_epi_moco,
    file_surf,
    delta_te=1.02,
    smooth=2.5,
    udir="y-",
    bw=16.304,
    nerode=1,
    cleanup=True,
):
    """This function computes a deformation field from a fieldmap acquisition and
    applies the inverse transformation to the undistorted surface. The following steps
    are performed:
        1. get median time series
        2. skullstrip epi
        3. register fieldmap to epi
        4. mask fieldmap
        5. prepare field
        6. get deforamtion field
        7. apply inverse deformation to surfaces.
        8. remove intermediate files (optional).

    To run the script, FSL and Freesurfer have to be in the PATH environment. The
    basenames of the surface files should be in freesurfer convention with the
    hemisphere indicated as prefix.

    Parameters
    ----------
    file_fmap_magn : str
        Fieldmap magnitude image.
    file_fmap_phase : str
        Fieldmap phase difference image.
    file_epi : str
        Filename of raw time series.
    file_epi_moco : str
        Filname of motion corrected time series.
    file_surf : list
        List of surface filnames.
    delta_te : float, optional
        Echo time difference of fieldmap in ms. The default is 1.02.
    smooth : float, optional
        Smoothing kernel for fieldmap unmasking. The default is 2.5.
    udir : str, optional
        Direction for fieldmap unmasking. The default is "y-".
    bw : float, optional
        BandwidthPerPixelPhaseEncode in Hz/px. The default is 16.304.
    nerode : int, optional
        Number of skullstrip mask eroding iterations. The default is 1.
    cleanup : bool, optional
        Removes temporary files at the end of the script. The default is True.
    """
    # prepare path and filename
    path_fmap0, name_fmap0, ext_fmap0 = get_filename(file_fmap_magn)
    path_fmap1, name_fmap1, ext_fmap1 = get_filename(file_fmap_phase)
    _, name_data, ext_data = get_filename(file_epi)
    path_udata, name_udata, ext_udata = get_filename(file_epi_moco)

    # filename with file extension
    name_fmap0 += ext_fmap0
    name_fmap1 += ext_fmap1
    name_data += ext_data
    name_udata += ext_udata

    # change directory to fieldmap directory
    os.chdir(path_fmap0)

    # get matrix size in phase encoding direction from uncorrected epi
    data = nb.load(file_epi)
    phase_encode = data.header.get_dim_info()[1]
    image_matrix_phase_encode = data.header["dim"][phase_encode + 1]

    # calculate median epi
    udata = nb.load(file_epi_moco)
    arr_udata = udata.get_fdata()
    arr_udata_median = np.median(arr_udata, axis=3)
    udata_median = nb.Nifti1Image(arr_udata_median, udata.affine, udata.header)
    udata_median.header["dim"][0] = 3
    udata_median.header["dim"][4] = 1
    nb.save(udata_median, os.path.join(path_udata, "median_" + name_udata))

    # calculate skullstrip mask of that image
    skullstrip_epi(
        os.path.join(path_udata, "median_" + name_udata),
        roi_size=10,
        scale=0.75,
        nerode=1,
        ndilate=2,
        savemask=True,
        cleanup=True,
    )

    # erode skullstrip mask
    for _ in range(nerode):
        erode_fsl(
            os.path.join(path_udata, "mask_median_" + name_udata),
            os.path.join(path_udata, "mask_median_" + name_udata),
        )

    # register fmap1 to median epi (FLIRT)
    flirt(
        file_fmap_magn,
        os.path.join(path_udata, "median_" + name_udata),
        os.path.join(path_fmap0, "r" + name_fmap0),
        os.path.join(path_fmap0, "fmap2epi.txt"),
        cost_func="mutualinfo",
        interp_method="trilinear",
    )

    # apply registration to fmap2
    apply_flirt(
        file_fmap_phase,
        os.path.join(path_udata, "median_" + name_udata),
        os.path.join(path_fmap0, "fmap2epi.txt"),
        os.path.join(path_fmap1, "r" + name_fmap1),
        "trilinear",
    )

    # apply skullstrip mask to fmap1 and fmap2 and save with same header information
    fmap1_img = nb.load(os.path.join(path_fmap0, "r" + name_fmap0))
    arr_fmap1 = fmap1_img.get_fdata()
    fmap2_img = nb.load(os.path.join(path_fmap1, "r" + name_fmap1))
    arr_fmap2 = fmap2_img.get_fdata()
    mask_img = nb.load(os.path.join(path_udata, "mask_median_" + name_udata))
    arr_mask = mask_img.get_fdata()

    arr_fmap1 = arr_fmap1 * arr_mask
    arr_fmap2 = arr_fmap2 * arr_mask
    arr_fmap2 = arr_fmap2 + np.abs(np.min(arr_fmap2))
    arr_fmap2 = (
        arr_fmap2 / np.max(arr_fmap2) * 4095
    )  # rescale phase image to be within 0-4095

    fmap1_img = nb.Nifti1Image(arr_fmap1, fmap1_img.affine, fmap1_img.header)
    nb.save(fmap1_img, os.path.join(path_fmap0, "pr" + name_fmap0))
    fmap2_img = nb.Nifti1Image(arr_fmap2, fmap1_img.affine, fmap1_img.header)
    nb.save(fmap2_img, os.path.join(path_fmap1, "pr" + name_fmap1))

    # prepare fieldmap (saves fieldmap in rad/s)
    prepare_fieldmap_fsl(
        os.path.join(path_fmap0, "pr" + name_fmap0),
        os.path.join(path_fmap1, "pr" + name_fmap1),
        os.path.join(path_fmap0, "fieldmap.nii"),
        delta_te,
    )

    # effective echo spacing in s
    dwell_time = 1 / (bw * image_matrix_phase_encode)

    # unmask fieldmap (fsl.FUGUE)
    fugue_fsl(
        os.path.join(path_udata, name_udata),
        os.path.join(path_fmap0, "fieldmap.nii"),
        os.path.join(path_fmap0, "vdm.nii"),
        dwell_time,
        smooth,
        udir,
    )

    # warp coordinate mapping
    generate_coordinate_mapping(
        file_epi, 0, path_fmap0, suffix="fmap", time=False, write_output=True
    )

    # apply inverse fieldmap to coordinate mapping
    apply_fugue(
        os.path.join(path_fmap0, "cmap_fmap.nii"),
        os.path.join(path_fmap0, "vdm.nii"),
        udir,
        False,
    )

    # apply cmap to surface
    for f in file_surf:
        path_surf, _, _ = get_filename(f)
        deform_surface(
            input_surf=f,
            input_orig=os.path.join(path_udata, "median_" + name_udata),
            input_deform=os.path.join(path_fmap0, "cmap_fmap_unwarped.nii"),
            input_target=os.path.join(path_udata, "median_" + name_udata),
            path_output=path_surf,
            input_mask=None,
            interp_method="trilinear",
            smooth_iter=0,
            flip_faces=False,
            cleanup=True,
        )

    # delete created files
    if cleanup:
        os.remove(os.path.join(path_fmap0, "cmap_fmap.nii"))
        os.remove(os.path.join(path_fmap0, "cmap_fmap_unwarped.nii"))
        os.remove(os.path.join(path_fmap0, "fieldmap.nii"))
        os.remove(os.path.join(path_fmap0, "fmap2epi.txt"))
        os.remove(
            os.path.join(path_fmap1, os.path.splitext(name_fmap1)[0] + "_flirt.mat")
        )
        os.remove(os.path.join(path_fmap0, "r" + name_fmap0))
        os.remove(os.path.join(path_fmap0, "pr" + name_fmap0))
        os.remove(os.path.join(path_fmap1, "r" + name_fmap1))
        os.remove(os.path.join(path_fmap1, "pr" + name_fmap1))
        os.remove(
            os.path.join(path_fmap0, os.path.splitext(name_udata)[0]) + "_unwarped.nii"
        )
        os.remove(os.path.join(path_fmap0, "vdm.nii"))
        os.remove(os.path.join(path_udata, "mask_median_" + name_udata))
        os.remove(os.path.join(path_udata, "median_" + name_udata))
        os.remove(os.path.join(path_udata, "pmedian_" + name_udata))
