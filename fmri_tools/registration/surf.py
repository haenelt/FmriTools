# -*- coding: utf-8 -*-
"""Deform surface mesh."""

import datetime
import os
import shutil as sh
import sys

import nibabel as nb
import numpy as np
from nibabel.freesurfer.io import read_geometry, write_geometry

from ..io.affine import read_vox2ras_tkr
from ..io.filename import get_filename
from ..io.vol import mri_convert
from ..preprocessing.fieldmap import fugue_fsl, prepare_fieldmap_fsl
from ..registration.cmap import generate_coordinate_mapping
from ..registration.fsl import flirt
from ..registration.mapping import mri_vol2surf
from ..registration.transform import apply_affine_chunked, apply_flirt, apply_fugue
from ..segmentation.skullstrip import skullstrip_epi
from ..surface.mesh import Mesh
from ..surface.smooth import mris_smooth
from ..utils.roi import erode_fsl

__all__ = ["apply_fieldmap", "deform_surface"]


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

    Returns
    -------
    None.

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

    Returns
    -------
    None.

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
