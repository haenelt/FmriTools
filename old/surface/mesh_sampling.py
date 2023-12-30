# -*- coding: utf-8 -*-

import datetime
import os
import shutil as sh
import sys

import nibabel as nb
import numpy as np
from sh import gunzip

from ..io.filename import get_filename
from ..io.surf import write_mgh
from ..mapping.map2surface import map2surface
from ..registration.cmap import generate_coordinate_mapping
from ..surface.deform_surface import deform_surface
from ..utils.resample_volume import resample_volume


def _rescale_cmap(file_cmap, dim, dim_upsampled):
    """Rescales and overwrites coordinate mapping if upsampled."""

    # load cmap
    cmap = nb.load(file_cmap)
    arr_cmap = cmap.get_fdata()

    # rescale
    for i in range(cmap.header["dim"][4]):
        arr_cmap[:, :, :, i] = arr_cmap[:, :, :, i] / dim[i] * dim_upsampled[i]

    # overwrite cmap
    cmap = nb.Nifti1Image(arr_cmap, cmap.affine, cmap.header)
    nb.save(cmap, file_cmap)


def mesh_sampling(
    surf_in,
    vol_in,
    write_output=False,
    path_output="",
    source2target_in="",
    interp_method="nearest",
    r=[0.4, 0.4, 0.4],
    interp_upsample="Cu",
    cleanup=True,
):
    """Mesh sampling.

    This function samples data onto a surface mesh. Optionally, the volume can
    be upsampled and a coordinate mapping can be applied to transform the
    surface mesh to the space of the input volume.

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
    r :list, optional
        Destination voxel size after upsampling (performed if not None). The
        default is [0.4,0.4,0.4].
    interp_upsample : str, optional
        Interpolation method if volume is upsampled. Possible arguments are NN,
        Li, Cu and Bk. The default is "Cu".
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

    # upsample volumes and rescale cmap
    if r:
        _, _, ext_vol = get_filename(file_vol)
        _, _, ext_s2t = get_filename(file_s2t)

        resample_volume(
            file_vol,
            os.path.join(path_tmp, name_vol + "_upsampled" + ext_vol),
            dxyz=r,
            rmode=interp_upsample,
        )

        resample_volume(
            file_s2t,
            os.path.join(path_tmp, name_s2t + "_upsampled" + ext_s2t),
            dxyz=r,
            rmode="Linear",
        )

        file_vol = os.path.join(path_tmp, name_vol + "_upsampled" + ext_vol)
        file_s2t = os.path.join(path_tmp, name_s2t + "_upsampled" + ext_s2t)

        _rescale_cmap(
            file_s2t,
            nb.load(vol_in).header["dim"][1:4] - 1,
            nb.load(file_vol).header["dim"][1:4] - 1,
        )

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
