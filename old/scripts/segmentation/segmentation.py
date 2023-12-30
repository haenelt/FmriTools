# -*- coding: utf-8 -*-
"""
Full cortex reconstruction

The purpose of the following script is to compute the segmentation for hires
MP2RAGE data and bring the resulting surface data of the cortical boundaries
into a regular grid representation. The script is divided into five parts.
Single steps of each part are stated below:

The script needs an installation of freesurfer.

Part 1 (Run recon-all pipeline)
    *flat image background denoising
    *bias field correction of the flat image
    *volume thresholding
    *recon1 without skullstripping
    *compute skullstrip mask (inv2)
    *recon23
Part 2 (Manual correction of white surface)
    *how to apply manual corrections is explained in a separate readme.
Part 3 (Manual correction of pial surface)
    *how to apply manual corrections is explained in a separate readme.
Part 4
    *inward shift of white surface to counteract the bias in mp2rage data
    *smooth pial and white surface (optional)
    *compute upsampled surface mesh
    *compute morphological data onto upsampled surface mesh
    *compute volumetric layers (optional)
Part 5
    *surface flattening of occipital patch
    *orthographic projection of flattened patch
    *visualise morphological data
    *visualise distortion
    *how to define a patch is explained in a separate readme.

"""

import datetime
import os

from cortex.polyutils import Surface
from nibabel.freesurfer.io import read_geometry
from nipype.interfaces.freesurfer import ApplyVolTransform

from ..layer.calc_equivol_surf import calc_equivol_surf
from ..mapping.map2grid import map2grid
from ..mapping.morph2dense import morph2dense
from ..matlab import MatlabCommand
from ..segmentation.get_ribbon_fsurf import get_ribbon_fsurf
from ..segmentation.get_thickness_fsurf import get_thickness_fsurf
from ..segmentation.include_pial_correction import include_pial_correction
from ..segmentation.orthographic_projection import orthographic_projection
from ..segmentation.shift_white import shift_white
from ..surface.flattening import surface_flattening
from ..surface.get_curvature import get_curvature
from ..surface.smooth import mris_smooth
from ..surface.upsample_surf_mesh import upsample_surf_mesh
from ..utils.bias import robust_combination
from ..utils.multiply_images import multiply_images
from ..utils.volume_threshold import volume_threshold

# input data
fileUNI = "/data/pt_01880/Experiment1_ODC/p1/anatomy/S6_MP2RAGE_0p7_UNI_Images_2.45.nii"
fileINV1 = "/data/pt_01880/Experiment1_ODC/p1/anatomy/S3_MP2RAGE_0p7_INV1_2.45.nii"
fileINV2 = "/data/pt_01880/Experiment1_ODC/p1/anatomy/S4_MP2RAGE_0p7_INV2_2.45.nii"
pathEXPERT = "/data/hu_haenelt/projects/FmriTools/scripts/segmentation"
namePATCH = "occip"
sub = "freesurfer"
part = 1

# parameters
reg_background = 8  # parameter for background noise removal (part 1)
w_shift = -0.5  # white surface shift (part 4)
niter_smooth = 2  # number of smoothing iterations for white and pial surface (part 4)
niter_upsample = 1  # number of upsampling iterations (part 4)
method_upsample = "linear"  # upsampling method (part 4)
nsurf_layer = 0  # number of equivolumetric layers, set 0 to omit (part 4)
factor_layer = 0  # smoothing of area surfaces (part 4)
niter_layer = 0  # number of smoothing iterations (part 4)
imres_ortho = 0.25  # isotropic image resolution of the regular grid in mm (part 5)
theta_ortho = [0, 0]  # rotation of the regular grid in deg for each hemisphere (part 5)
alpha_ortho = 2  # alpha shape value for concave hull computation (part 5)
buffer_ortho = 0  # smooth out of concave hull (part 5)
sigma_map = 0.5  # isotropic smoothing of data onto regular grid (part 5)

# do not edit below

# parameters
hemi = ["lh", "rh"]

# get path and file name
path_split = os.path.split(fileUNI)
path = path_split[0]
file = path_split[1]

# set folder structure
path_bias = os.path.join(path, "bias")
path_skull = os.path.join(path, "skull")
path_dense = os.path.join(path, "dense")
path_layer = os.path.join(path, "layer")
path_ortho = os.path.join(path, "ortho")

# write log
fileID = open(os.path.join(path, "segmentation_info.txt"), "a")
fileID.write(
    "Part "
    + str(part)
    + ": "
    + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    + "\n"
)
fileID.close()

if part == 1:
    # background noise removal
    print("Background noise removal")
    robust_combination(fileUNI, fileINV1, fileINV2, reg_background, path_bias)

    # bias field correction
    print("Bias field correction")
    matlab = MatlabCommand(
        "ft_bias_field_correction", os.path.join(path_bias, "n" + file)
    )
    matlab.run()

    # volume threshold
    print("Volume threshold")
    volume_threshold(os.path.join(path_bias, "mn" + file), "", 4095)

    # autorecon1 without skullstrip removal
    print("Autorecon1")
    os.system(
        "recon-all"
        + " -i "
        + os.path.join(path_bias, "mn" + file)
        + " -hires"
        + " -autorecon1"
        + " -noskullstrip"
        + " -sd "
        + path
        + " -s "
        + sub
        + " -parallel"
    )

    # skullstrip anatomy
    print("Skullstrip INV2")
    matlab = MatlabCommand("ft_skullstrip_spm12", fileINV2, path)
    matlab.run()

    # bring skullstrip_mask in conformed space (mri_vol2vol, NN)
    transmask = ApplyVolTransform()
    transmask.inputs.source_file = os.path.join(path_skull, "skullstrip_mask.nii")
    transmask.inputs.target_file = os.path.join(path, sub, "mri", "orig.mgz")
    transmask.inputs.reg_header = True
    transmask.inputs.interp = "nearest"
    transmask.inputs.transformed_file = os.path.join(
        path, sub, "mri", "skullstrip_mask.nii"
    )
    transmask.inputs.args = "--no-save-reg"
    transmask.run()

    # apply skullstrip mask to T1 and save as brainmask
    multiply_images(
        os.path.join(path, sub, "mri", "T1.mgz"),
        os.path.join(path, sub, "mri", "skullstrip_mask.nii"),
        os.path.join(path, sub, "mri", "brainmask.mgz"),
    )

    multiply_images(
        os.path.join(path, sub, "mri", "T1.mgz"),
        os.path.join(path, sub, "mri", "skullstrip_mask.nii"),
        os.path.join(path, sub, "mri", "brainmask.auto.mgz"),
    )

    # autorecon2 and autorecon3
    print("Autorecon2 and Autorecon3")
    os.system(
        "recon-all"
        + " -hires"
        + " -autorecon2"
        + " -autorecon3"
        + " -sd "
        + path
        + " -s "
        + sub
        + " -expert "
        + os.path.join(pathEXPERT, "expert.opts")
        + " -xopts-overwrite"
        + " -parallel"
    )

    # write log
    fileID = open(os.path.join(path, "segmentation_info.txt"), "a")
    fileID.write("Regularisation for background removal: " + str(reg_background) + "\n")
    fileID.close()

elif part == 2:
    # GM/WM surface correction using modified wm.mgz
    print("Autorecon2-wm")
    os.system(
        "recon-all"
        + " -hires"
        + " -autorecon2-wm"
        + " -autorecon3"
        + " -sd "
        + path
        + " -s "
        + sub
        + " -expert "
        + os.path.join(pathEXPERT, "expert.opts")
        + " -xopts-overwrite"
        + " -parallel"
    )

elif part == 3:
    # Convert manual corrections in pial_edit.mgz to brain.finalsurfs.manedit.mgz
    print("Create brain.finalsurfs.manedit.mgz")
    include_pial_correction(path, sub)

    # GM/CSF surface correction using modified brain.finalsurfs.manedit.mgz
    print("Autorecon-pial")
    os.system(
        "recon-all"
        + " -hires"
        + " -autorecon-pial"
        + " -sd "
        + path
        + " -s "
        + sub
        + " -expert "
        + os.path.join(pathEXPERT, "expert.opts")
        + " -xopts-overwrite"
        + " -parallel"
    )

elif part == 4:
    # output folder (freesurfer trash folder)
    path_trash = os.path.join(path, sub, "trash")

    # get date string for moved files
    date = datetime.datetime.now().strftime("%Y%m%d%H%M")

    # inward shift of final white surface
    print("Finalise white surface")
    shift_white(path, sub, w_shift)

    if niter_smooth != 0:
        print("Smooth white and pial surface")
        file_surf = ["white", "pial"]
        for i in range(len(file_surf)):
            for j in range(len(hemi)):
                file_in = os.path.join(path, sub, "surf", hemi[j] + "." + file_surf[i])
                file_out = os.path.join(
                    path_trash, hemi[j] + "." + file_surf[i] + "_backup_" + date
                )
                os.rename(file_in, file_out)
                mris_smooth(file_out, file_in, niter_smooth)

    # generate new curvature, thickness and ribbon files
    print("Compute new morphological files")
    for i in range(len(hemi)):
        file_in = os.path.join(path, sub, "surf", hemi[i] + ".curv")
        file_out = os.path.join(path_trash, hemi[i] + ".curv_backup_" + date)
        os.rename(file_in, file_out)
        get_curvature(
            os.path.join(path, sub, "surf", hemi[i] + ".white"),
            os.path.join(path, sub, "surf"),
        )

    get_thickness_fsurf(path, sub)
    get_ribbon_fsurf(path, sub)

    # upsample surface mesh
    print("Upsample surface mesh")
    orig_params = []
    dense_params = []
    file_surf = ["sphere", "white", "pial", "inflated"]  # list of surfaces to subdivide
    for i in range(len(file_surf)):
        for j in range(len(hemi)):
            file_in = os.path.join(path, sub, "surf", hemi[j] + "." + file_surf[i])
            file_out = os.path.join(path_dense, hemi[j] + "." + file_surf[i])
            upsample_surf_mesh(file_in, file_out, niter_upsample, method_upsample)

            if i == 0:
                vtx, fac = read_geometry(file_in)
                vtx_dense, fac_dense = read_geometry(file_out)
                orig = [len(vtx[:, 0]), Surface(vtx, fac).avg_edge_length]
                dense = [
                    len(vtx_dense[:, 0]),
                    Surface(vtx_dense, fac_dense).avg_edge_length,
                ]
                orig_params.extend(orig)
                dense_params.extend(dense)

    # transform curv to dense surface
    print("Transform morphological files to dense")
    for i in range(len(hemi)):
        morph2dense(
            os.path.join(path, sub, "surf", hemi[i] + ".sphere"),
            os.path.join(path_dense, hemi[i] + ".sphere"),
            os.path.join(path, sub, "surf", hemi[i] + ".curv"),
            path_dense,
        )

        morph2dense(
            os.path.join(path, sub, "surf", hemi[i] + ".sphere"),
            os.path.join(path_dense, hemi[i] + ".sphere"),
            os.path.join(path, sub, "surf", hemi[i] + ".thickness"),
            path_dense,
        )

    # calculate volumetric surfaces
    if nsurf_layer != 0:
        print("Compute volumetric layers")
        for i in range(len(hemi)):
            file_white = os.path.join(path_dense, hemi[i] + ".white")
            file_pial = os.path.join(path_dense, hemi[i] + ".pial")
            calc_equivol_surf(
                file_white,
                file_pial,
                nsurf_layer,
                factor_layer,
                niter_layer,
                hemi[i],
                path_layer,
            )

    # write log
    fileID = open(os.path.join(path, "segmentation_info.txt"), "a")
    fileID.write("Inward shift of white surface: " + str(w_shift) + "\n")
    fileID.write("Number of smoothing iterations: " + str(niter_smooth) + "\n")
    fileID.write("Number of upsampling iterations: " + str(niter_upsample) + "\n")
    fileID.write("Upsampling method: " + method_upsample + "\n")
    fileID.write(
        "Number of nodes in original surface (left): " + str(orig_params[0]) + "\n"
    )
    fileID.write(
        "Average edge length in original surface (left): " + str(orig_params[1]) + "\n"
    )
    fileID.write(
        "Number of nodes in original surface (right): " + str(orig_params[2]) + "\n"
    )
    fileID.write(
        "Average edge length in original surface (right): " + str(orig_params[3]) + "\n"
    )
    fileID.write(
        "Number of nodes in dense surface (left): " + str(dense_params[0]) + "\n"
    )
    fileID.write(
        "Average edge length in dense surface (left): " + str(dense_params[1]) + "\n"
    )
    fileID.write(
        "Number of nodes in dense surface (right): " + str(dense_params[2]) + "\n"
    )
    fileID.write(
        "Average edge length in dense surface (right): " + str(dense_params[3]) + "\n"
    )
    fileID.write("Number of volumetric surfaces: " + str(nsurf_layer) + "\n")
    fileID.write("Smoothing factor for layering: " + str(factor_layer) + "\n")
    fileID.write(
        "Number of smoothing iterations for layering: " + str(niter_layer) + "\n"
    )
    fileID.close()

elif part == 5:
    # surface flattening
    print("Surface flattening")
    for i in range(len(hemi)):
        surface_flattening(
            os.path.join(path_dense, hemi[i] + ".white"),
            os.path.join(path_dense, hemi[i] + "." + namePATCH + ".patch.3d"),
            path_dense,
            cleanup=True,
        )

    # regular grid interpolation
    print("Orthographic projection")
    nvoxel_params = []
    ind_ratio_params = []
    for i in range(len(hemi)):
        nvoxel, ind_ratio = orthographic_projection(
            os.path.join(path_dense, hemi[i] + "." + namePATCH + ".patch.flat"),
            imres_ortho,
            theta_ortho[i],
            alpha_ortho,
            buffer_ortho,
            path_ortho,
        )
        nvoxel_params.append(nvoxel)
        ind_ratio_params.append(ind_ratio)

    # map distortion data onto grid
    print("Map distortion data onto grid")
    for i in range(len(hemi)):
        map2grid(
            os.path.join(
                path_ortho, hemi[i] + "." + namePATCH + ".patch.flat.cmap.nii"
            ),
            os.path.join(path_dense, hemi[i] + ".curv"),
            sigma_map,
            path_ortho,
            hemi[i] + "." + namePATCH + ".curv",
        )
        map2grid(
            os.path.join(
                path_ortho, hemi[i] + "." + namePATCH + ".patch.flat.cmap.nii"
            ),
            os.path.join(path_dense, hemi[i] + ".thickness"),
            sigma_map,
            path_ortho,
            hemi[i] + "." + namePATCH + ".thickness",
        )

    # write log
    fileID = open(os.path.join(path, "segmentation_info.txt"), "a")
    fileID.write(namePATCH + ": Image resolution of grid -> " + str(imres_ortho) + "\n")
    fileID.write(namePATCH + ": Grid rotation (left) -> " + str(theta_ortho[0]) + "\n")
    fileID.write(namePATCH + ": Grid rotation (right) -> " + str(theta_ortho[1]) + "\n")
    fileID.write(
        namePATCH + ": Alpha shape for concave hull -> " + str(alpha_ortho) + "\n"
    )
    fileID.write(namePATCH + ": Concave hull buffer -> " + str(buffer_ortho) + "\n")
    fileID.write(
        namePATCH
        + ": Number of grid voxels within patch (left) -> "
        + str(nvoxel_params[0])
        + "\n"
    )
    fileID.write(
        namePATCH
        + ": Ratio of unique grid indices (left) -> "
        + str(ind_ratio_params[0])
        + "\n"
    )
    fileID.write(
        namePATCH
        + ": Number of grid voxels within patch (right) -> "
        + str(nvoxel_params[1])
        + "\n"
    )
    fileID.write(
        namePATCH
        + ": Ratio of unique grid indices (right) -> "
        + str(ind_ratio_params[1])
        + "\n"
    )
    fileID.write(namePATCH + ": Grid mapping (sigma) -> " + str(sigma_map) + "\n")
    fileID.close()

else:
    print("Which part do you want to run?")
