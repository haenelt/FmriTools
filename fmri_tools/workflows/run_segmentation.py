# -*- coding: utf-8 -*-
"""Full cortex reconstruction

The purpose of the following workflow is to compute the segmentation for hires MP2RAGE
(or MEMPRAGE) data and bring the resulting surface data of the cortical boundaries into
a regular grid representation. The script is divided into five parts. Single steps of of
each part are stated below:

The workflow needs an installation of freesurfer.

Part 1 (Run recon-all pipeline)
    *flat image background denoising (only for MP2RAGE)
    *bias field correction of the flat image
    *volume thresholding
    *recon1 without skullstripping
    *compute skullstrip mask (MP2RAGE: inv2, MEMPRAGE: T1)
    *recon23
Part 2 (Manual correction of white surface)
    *how to apply manual corrections is explained further below.
Part 3 (Manual correction of pial surface)
    *how to apply manual corrections is explained further below.
Part 4
    *inward shift of white surface to counteract the bias in mp2rage data (not applied)
    *smooth pial and white surface (optional)
    *compute upsampled surface mesh
    *compute morphological data onto upsampled surface mesh
    *compute volumetric layers (optional)
Part 5
    *surface flattening of occipital patch
    *orthographic projection of flattened patch
    *visualise morphological data
    *visualise distortion
    *how to define a patch is explained further below.

HOWTO: manual white surface correction
    *copy old `wm.mgz` in freesurfer trash folder and name like
    `wm_backup_201812061015.mgz`
    *copy old white surfaces the same way
    *edit `wm.mgz` (brush value: 255, eraser value: 1)
    *best to overlay with orig and surfaces

HOWTO: manual pial surface correction
    *move old `brain.finalsurfs.mgz` to trash folder as in part 2 (remove from `mri`
    folder!)
    *copy old pial and white surface the same way
    *copy `orig.mgz` and save as `pial_edit.mgz` in the same folder
    *`mri_convert pial_edit.mgz pial_edit.mgz -odt float`
    *apply changes (brush value: 256, eraser value: -1)
    *best to overlay with surfaces
    *N.B. in freesurfer, manual changes of the pial surfaces are applied by creating the
    file `brain.finalsurfs.manedit.mgz` (copy of `brainmask.mgz`, brush value: 255,
    eraser: 1). This file will be created when the script is run.

HOWTO: defining a patch for surface flattening
    *open `tksurfer` and define manually the patch of the occipital pole
    *if you define the patch onto the upsampled surface mesh, you have to open first the
    original surface in `tksurfer` and then load the upsampled surface in the opened GUI
    *load the inflated surface
    *rotate to the medial surface
    *select points along the calcarine fissure and press the button "Cut line"
    *select 3 points to define the cutting plane: 2 on medial side and 1 on lateral side
    *choose a 4th points to specify which portion of surface to keep and press button
    "Cut plane"
    *save file: File > Patch > Save as file `<HEMI>.<name_patch>.patch.3d`
    *save the file in the dense SUBfolder
    *after flattening you can visualise the patch by loading first the inflated surface
    in tksurfer, then File > Patch > Load patch ...
"""

import datetime
import os
import subprocess
from argparse import SUPPRESS, ArgumentParser

import fmri_tools
from nibabel.freesurfer.io import read_geometry

from ..io.filename import get_filename
from ..matlab import MatlabCommand
from ..registration.mapping import map2grid, morph2dense
from ..registration.transform import apply_header
from ..segmentation.flat import orthographic_projection, surface_flattening
from ..segmentation.layer import calc_equivol_surf
from ..segmentation.mask import mask_ribbon
from ..segmentation.surf import (
    mris_curvature,
    mris_thickness,
    shift_white,
    upsample_surf_mesh,
)
from ..segmentation.vol import include_pial_correction, robust_combination
from ..surface.mesh import Mesh
from ..surface.smooth import mris_smooth
from ..utils.calc import multiply_images, volume_threshold

# parameters
SUB = "freesurfer"
HEMI = ["lh", "rh"]
PATH_EXPERT = os.path.join(os.path.dirname(fmri_tools.__file__), "segmentation")
REG_BACKGROUND = 8  # parameter for background noise removal (part 1)
W_SHIFT = 0.0  # white surface shift (part 4)
NITER_SMOOTH = 2  # number of smoothing iterations for white and pial surface (part 4)
NITER_UPSAMPLE = 1  # number of upsampling iterations (part 4)
METHOD_UPSAMPLE = "linear"  # upsampling method (part 4)
IMRES_ORTHO = 0.25  # isotropic image resolution of the regular grid in mm (part 5)
THETA_ORTHO = [0, 0]  # rotation of the regular grid in deg for each hemisphere (part 5)
ALPHA_ORTHO = 2  # alpha shape value for concave hull computation (part 5)
BUFFER_ORTHO = 0  # smooth out of concave hull (part 5)
SIGMA_MAP = 0.5  # isotropic smoothing of data onto regular grid (part 5)


def _get_parser():
    """Parse command line inputs.

    Returns
    -------
    parser.parse_args() : argparse dict
    """
    # Disable default help
    parser = ArgumentParser(add_help=False)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    # Add back help
    optional.add_argument(
        "--help",
        action="help",
        default=SUPPRESS,
        help="Show this help message and exit.",
    )
    required.add_argument(
        "--uni",
        dest="uni",
        type=str,
        help="File name of UNI of T1 image.",
        required=True,
    )
    required.add_argument(
        "--part",
        dest="part",
        type=int,
        help="Which part to run.",
        required=True,
    )
    optional.add_argument(
        "--inv1",
        dest="inv1",
        type=str,
        help="File name of INV1 image",
        default="",
    )
    optional.add_argument(
        "--inv2",
        dest="inv2",
        type=str,
        help="File name of INV2 image",
        default="",
    )
    optional.add_argument(
        "--flair",
        dest="flair",
        type=str,
        help="File name of FLAIR image (default: %(default)s)",
        default="",
    )
    optional.add_argument(
        "--patch",
        dest="name_patch",
        type=str,
        help="Name of flat patch (default: %(default)s)",
        default="occip",
    )
    optional.add_argument(
        "--layer",
        dest="n_layer",
        type=int,
        help=(
            "Number of equivolumetric layers, set 0 to omit. This parameter is only \
                relevant in part 4 (default: %(default)s)."
        ),
        default=11,
    )

    return parser


def segmentation_workflow(uni, part, inv1, inv2, flair, name_patch, n_layer):
    """Workflow for MP2RAGE segmentation.

    Parameters
    ----------
    uni : str
        File name of MP2RAGE UNI or MEMPRAGE T1 image.
    part : int
        Part of segmentation that should be run.
    inv1 : str, optional
        File name of MP2RAGE INV1 image.
    inv2 : str, optional
        File name of MP2RAGE INV2 image.
    flair : str, optional
        File name of FLAIR image.
    name_patch : str, optional
        Name of flat patch.
    n_layer : int, optional
        Number of layers.

    Raises
    ------
    ValueError
        If incorrect part is entered.
    """
    # check part
    if part < 1 or part > 5:
        raise ValueError("part must be an integer between 1 and 5.")

    # get path and file name
    path_uni, name_uni, ext_uni = get_filename(uni)
    file_uni = f"{name_uni}{ext_uni}"

    # set folder structure
    path_bias = os.path.join(path_uni, "bias")
    path_skull = os.path.join(path_uni, "skull")
    path_dense = os.path.join(path_uni, "dense")
    path_layer = os.path.join(path_uni, "layer")
    path_ortho = os.path.join(path_uni, "ortho")

    # write log
    time_stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    file_id = open(
        os.path.join(path_uni, "segmentation_info.txt"), "a", encoding="utf-8"
    )
    file_id.write(f"Part {part}: {time_stamp}\n")
    file_id.close()

    if part == 1:
        # background noise removal
        if not os.path.exists(path_bias):
            os.mkdir(path_bias)
        if inv1 and inv2:
            print("Background noise removal")
            robust_combination(uni, inv1, inv2, REG_BACKGROUND, path_bias)
        else:
            os.symlink(uni, os.path.join(path_bias, "n" + file_uni))

        # bias field correction
        print("Bias field correction")
        matlab = MatlabCommand(
            "ft_bias_field_correction", os.path.join(path_bias, "n" + file_uni)
        )
        matlab.run()

        # volume threshold
        print("Volume threshold")
        volume_threshold(os.path.join(path_bias, "mn" + file_uni), "", 4095)

        # autorecon1 without skullstrip removal
        print("Autorecon1")
        command = "recon-all"
        command += f" -i {os.path.join(path_bias, 'mn' + file_uni)}"
        command += " -hires -autorecon1 -noskullstrip"
        command += f" -sd {path_uni}"
        command += f" -s {SUB}"
        command += " -parallel"
        if flair:
            command += f" -FLAIR {flair} -FLAIRpial"

        print("Execute: " + command)
        try:
            subprocess.run([command], shell=True, check=False)
        except subprocess.CalledProcessError:
            print("Execuation failed!")

        # skullstrip anatomy
        if inv2:
            print("Skullstrip INV2")
            matlab = MatlabCommand("ft_skullstrip_spm12", inv2, path_uni)
        else:
            print("Skullstrip T1")
            matlab = MatlabCommand("ft_skullstrip_spm12", uni, path_uni)
        matlab.run()

        # bring skullstrip_mask in conformed space (mri_vol2vol, NN)
        apply_header(
            os.path.join(path_skull, "skullstrip_mask.nii"),
            os.path.join(path_uni, SUB, "mri", "orig.mgz"),
            os.path.join(path_uni, SUB, "mri", "skullstrip_mask.nii"),
            interp_method="nearest",
        )

        # apply skullstrip mask to T1 and save as brainmask
        multiply_images(
            os.path.join(path_uni, SUB, "mri", "T1.mgz"),
            os.path.join(path_uni, SUB, "mri", "skullstrip_mask.nii"),
            os.path.join(path_uni, SUB, "mri", "brainmask.mgz"),
        )

        multiply_images(
            os.path.join(path_uni, SUB, "mri", "T1.mgz"),
            os.path.join(path_uni, SUB, "mri", "skullstrip_mask.nii"),
            os.path.join(path_uni, SUB, "mri", "brainmask.auto.mgz"),
        )

        # autorecon2 and autorecon3
        print("Autorecon2 and Autorecon3")
        command = "recon-all"
        command += " -hires -autorecon2 -autorecon3"
        command += f" -sd {path_uni}"
        command += f" -s {SUB}"
        command += f" -expert {os.path.join(PATH_EXPERT, 'expert.opts')}"
        command += " -xopts-overwrite -parallel"
        if flair:
            command += f" -FLAIR {flair} -FLAIRpial"

        print("Execute: " + command)
        try:
            subprocess.run([command], shell=True, check=False)
        except subprocess.CalledProcessError:
            print("Execuation failed!")

        # write log
        file_id = open(
            os.path.join(path_uni, "segmentation_info.txt"), "a", encoding="utf-8"
        )
        file_id.write(
            "Regularisation for background removal: " + str(REG_BACKGROUND) + "\n"
        )
        file_id.close()

    elif part == 2:
        # GM/WM surface correction using modified wm.mgz
        print("Autorecon2-wm")
        command = "recon-all"
        command += " -hires -autorecon2-wm -autorecon3"
        command += f" -sd {path_uni}"
        command += f" -s {SUB}"
        command += f" -expert {os.path.join(PATH_EXPERT, 'expert.opts')}"
        command += " -xopts-overwrite -parallel"
        if flair:
            command += f" -FLAIR {flair} -FLAIRpial"

        print("Execute: " + command)
        try:
            subprocess.run([command], shell=True, check=False)
        except subprocess.CalledProcessError:
            print("Execuation failed!")

    elif part == 3:
        # Convert manual corrections in pial_edit.mgz to brain.finalsurfs.manedit.mgz
        print("Create brain.finalsurfs.manedit.mgz")
        include_pial_correction(path_uni, SUB)

        # GM/CSF surface correction using modified brain.finalsurfs.manedit.mgz
        print("Autorecon-pial")
        command = "recon-all"
        command += " -hires -autorecon-pial"
        command += f" -sd {path_uni}"
        command += f" -s {SUB}"
        command += f" -expert {os.path.join(PATH_EXPERT, 'expert.opts')}"
        command += " -xopts-overwrite -parallel"
        if flair:
            command += f" -FLAIR {flair} -FLAIRpial"

        print("Execute: " + command)
        try:
            subprocess.run([command], shell=True, check=False)
        except subprocess.CalledProcessError:
            print("Execuation failed!")

    elif part == 4:
        # output folder (freesurfer trash folder)
        path_trash = os.path.join(path_uni, SUB, "trash")

        # get date string for moved files
        date = datetime.datetime.now().strftime("%Y%m%d%H%M")

        # inward shift of final white surface
        print("Finalise white surface")
        shift_white(path_uni, SUB, W_SHIFT)

        if NITER_SMOOTH != 0:
            print("Smooth white and pial surface")
            file_surf = ["white", "pial"]
            for _file in file_surf:
                for _h in HEMI:
                    file_in = os.path.join(path_uni, SUB, "surf", f"{_h}.{_file}")
                    file_out = os.path.join(path_trash, f"{_h}.{_file}_backup_{date}")
                    os.rename(file_in, file_out)
                    mris_smooth(file_out, file_in, NITER_SMOOTH)

        # generate new curvature, thickness and ribbon files
        print("Compute new morphological files")
        for _h in HEMI:
            file_in = os.path.join(path_uni, SUB, "surf", f"{_h}.curv")
            file_out = os.path.join(path_trash, f"{_h}.curv_backup_{date}")
            os.rename(file_in, file_out)
            mris_curvature(
                os.path.join(path_uni, SUB, "surf", f"{_h}.white"),
                os.path.join(path_uni, SUB, "surf"),
            )

        mris_thickness(path_uni, SUB)
        mask_ribbon(path_uni, SUB)

        # upsample surface mesh
        print("Upsample surface mesh")
        orig_params = []
        dense_params = []
        file_surf = [
            "sphere",
            "white",
            "pial",
            "inflated",
        ]  # list of surfaces to SUBdivide
        for i, _file in enumerate(file_surf):
            for _h in HEMI:
                file_in = os.path.join(path_uni, SUB, "surf", f"{_h}.{_file}")
                file_out = os.path.join(path_dense, f"{_h}.{_file}")
                upsample_surf_mesh(file_in, file_out, NITER_UPSAMPLE, METHOD_UPSAMPLE)

                if i == 0:
                    vtx, fac = read_geometry(file_in)
                    vtx_dense, fac_dense = read_geometry(file_out)
                    orig = [len(vtx[:, 0]), Mesh(vtx, fac).avg_edge_length]
                    dense = [
                        len(vtx_dense[:, 0]),
                        Mesh(vtx_dense, fac_dense).avg_edge_length,
                    ]
                    orig_params.extend(orig)
                    dense_params.extend(dense)

        # transform curv to dense surface
        print("Transform morphological files to dense")
        for _h in HEMI:
            morph2dense(
                os.path.join(path_uni, SUB, "surf", f"{_h}.sphere"),
                os.path.join(path_dense, f"{_h}.sphere"),
                os.path.join(path_uni, SUB, "surf", f"{_h}.curv"),
                path_dense,
            )

            morph2dense(
                os.path.join(path_uni, SUB, "surf", f"{_h}.sphere"),
                os.path.join(path_dense, f"{_h}.sphere"),
                os.path.join(path_uni, SUB, "surf", f"{_h}.thickness"),
                path_dense,
            )

        # calculate volumetric surfaces
        if n_layer != 0:
            print("Compute volumetric layers")
            for _h in HEMI:
                file_white = os.path.join(path_dense, f"{_h}.white")
                file_pial = os.path.join(path_dense, f"{_h}.pial")
                calc_equivol_surf(
                    file_white,
                    file_pial,
                    n_layer,
                    HEMI[i],
                    path_layer,
                )

        # write log
        file_id = open(
            os.path.join(path_uni, "segmentation_info.txt"), "a", encoding="utf-8"
        )
        file_id.write(f"Inward shift of white surface: {W_SHIFT}\n")
        file_id.write(f"Number of smoothing iterations: {NITER_SMOOTH}\n")

        file_id.write(f"Number of upsampling iterations: {NITER_UPSAMPLE}\n")
        file_id.write(f"Upsampling method: {METHOD_UPSAMPLE}\n")
        file_id.write(f"Number of nodes in original surface (left): {orig_params[0]}\n")
        file_id.write(
            f"Average edge length in original surface (left): {orig_params[1]}\n"
        )
        file_id.write(
            f"Number of nodes in original surface (right): {orig_params[2]}\n"
        )
        file_id.write(
            f"Average edge length in original surface (right): {orig_params[3]}\n"
        )
        file_id.write(f"Number of nodes in dense surface (left): {dense_params[0]}\n")
        file_id.write(
            f"Average edge length in dense surface (left): {dense_params[1]}\n"
        )
        file_id.write(f"Number of nodes in dense surface (right): {dense_params[2]}\n")
        file_id.write(
            f"Average edge length in dense surface (right): {dense_params[3]}\n"
        )
        file_id.write(f"Number of volumetric surfaces: {n_layer}\n")
        file_id.close()

    elif part == 5:
        # surface flattening
        print("Surface flattening")
        for _h in HEMI:
            surface_flattening(
                os.path.join(path_dense, f"{_h}.white"),
                os.path.join(path_dense, f"{_h}.{name_patch}.patch.3d"),
                path_dense,
                cleanup=True,
            )

        # regular grid interpolation
        print("Orthographic projection")
        nvoxel_params = []
        ind_ratio_params = []
        for _h in HEMI:
            nvoxel, ind_ratio = orthographic_projection(
                os.path.join(path_dense, f"{_h}.{name_patch}.patch.flat"),
                IMRES_ORTHO,
                THETA_ORTHO[i],
                ALPHA_ORTHO,
                BUFFER_ORTHO,
                path_ortho,
            )
            nvoxel_params.append(nvoxel)
            ind_ratio_params.append(ind_ratio)

        # map distortion data onto grid
        print("Map distortion data onto grid")
        for _h in HEMI:
            map2grid(
                os.path.join(path_ortho, f"{_h}.{name_patch}.patch.flat.cmap.nii"),
                os.path.join(path_dense, f"{_h}.curv"),
                SIGMA_MAP,
                path_ortho,
                f"{_h}.{name_patch}.curv",
            )
            map2grid(
                os.path.join(path_ortho, f"{_h}.{name_patch}.patch.flat.cmap.nii"),
                os.path.join(path_dense, f"{_h}.thickness"),
                SIGMA_MAP,
                path_ortho,
                f"{_h}.{name_patch}.thickness",
            )

        # write log
        file_id = open(
            os.path.join(path_uni, "segmentation_info.txt"), "a", encoding="utf-8"
        )
        file_id.write(f"{name_patch}: Image resolution of grid -> {IMRES_ORTHO}\n")
        file_id.write(f"{name_patch}: Grid rotation (left) -> {THETA_ORTHO[0]}\n")
        file_id.write(f"{name_patch}: Grid rotation (right) -> {THETA_ORTHO[1]}\n")
        file_id.write(f"{name_patch}: Alpha shape for concave hull -> {ALPHA_ORTHO}\n")
        file_id.write(f"{name_patch}: Concave hull buffer -> {BUFFER_ORTHO}\n")
        file_id.write(
            f"{name_patch}: Number of grid voxels within patch (left) -> {nvoxel_params[0]}\n"
        )
        file_id.write(
            f"{name_patch}: Ratio of unique grid indices (left) -> {ind_ratio_params[0]}\n"
        )
        file_id.write(
            f"{name_patch}: Number of grid voxels within patch (right) -> {nvoxel_params[1]}\n"
        )
        file_id.write(
            f"{name_patch}: Ratio of unique grid indices (right) -> {ind_ratio_params[1]}\n"
        )
        file_id.write(f"{name_patch}: Grid mapping (sigma) -> {SIGMA_MAP}\n")
        file_id.close()

    else:
        print("Which part do you want to run?")


def _main(argv=None):
    """MP2RAGE segmentation workflow."""
    options = _get_parser().parse_args(argv)
    kwargs = vars(options)
    segmentation_workflow(**kwargs)


if __name__ == "__main__":
    _main()
