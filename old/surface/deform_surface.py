# -*- coding: utf-8 -*-

import datetime
import os
import shutil as sh
import sys

import nibabel as nb
import numpy as np
from gbb.utils import remove_vertex
from nibabel.freesurfer.io import read_geometry, write_geometry

from ..io.affine import read_vox2ras_tkr
from ..io.filename import get_filename
from ..io.vol import mri_convert
from ..registration.mapping import mri_vol2surf
from ..registration.transform import apply_affine_chunked
from ..surface.smooth import mris_smooth


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
    """Deform surface using freesurfer.

    This function deforms a surface mesh in freesurfer convention using a
    coordinate map containing voxel coordinates. The computation takes quite a
    while because in the case of removed vertices, i.e. if a mask is given as
    input, the remaining faces are reindexed.

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
    _, name_orig, ext_orig = get_filename(input_orig)
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
    for i in range(len(components)):
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

        vtx_new, fac_new, ind_keep = remove_vertex(vtx_new, fac, ind_keep)

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
