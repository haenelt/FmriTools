# -*- coding: utf-8 -*-
"""ANTS registration."""

import math
import os
import shutil as sh
from glob import glob
from os.path import exists, join

import nibabel as nb
import numpy as np
from nibabel.freesurfer.io import write_geometry
from scipy.io import loadmat, savemat

from .. import check_installation, execute_command
from ..io.affine import apply_affine_chunked, read_vox2ras_tkr
from ..io.filename import get_filename
from ..io.vol import load_volume, save_volume
from ..matlab import MatlabCommand
from .cmap import generate_coordinate_mapping

# convenience labels
X = 0
Y = 1
Z = 2

__all__ = ["embedded_antsreg", "recursive_bbr"]


def embedded_antsreg(
    file_source,
    file_target,
    output_dir,
    run_rigid=True,
    rigid_iterations=1000,
    run_affine=False,
    affine_iterations=1000,
    run_syn=True,
    coarse_iterations=40,
    medium_iterations=50,
    fine_iterations=40,
    cost_function="MutualInformation",
    interpolation="NearestNeighbor",
    convergence=1e-6,
    ignore_affine=False,
    ignore_header=False,
):
    """Runs the rigid and/or Symmetric Normalization (SyN) algorithm of ANTs and formats
    the output deformations into voxel coordinate mappings as used in CBSTools
    registration and transformation routines. Three files are writte to disk:
    (1) transformed source file (suffix: <file_source>_ants-def), (2) forward cmap
    (suffix: <file_source>_ants-map) and (3) inverse cmap (suffix:
    <file_source>_ants-invmap) with file extension of the source image.

    Parameters
    ----------
    file_source: str
        File name or niimg of source image.
    file_target: str
        File name or niimg of reference image.
    output_dir: str
        Path to desired output directory, will be created if it does not exist.
    run_rigid: bool
        Whether or not to run a rigid registration first (default is False).
    rigid_iterations: float
        Number of iterations in the rigid step (default is 1000).
    run_affine: bool
        Whether or not to run a affine registration first (default is False).
    affine_iterations: float
        Number of iterations in the affine step (default is 1000).
    run_syn: bool
        Whether or not to run a SyN registration (default is True).
    coarse_iterations: float
        Number of iterations at the coarse level (default is 40).
    medium_iterations: float
        Number of iterations at the medium level (default is 50).
    fine_iterations: float
        Number of iterations at the fine level (default is 40).
    cost_function: {'CrossCorrelation', 'MutualInformation'}
        Cost function for the registration (default is 'MutualInformation').
    interpolation: {'NearestNeighbor', 'Linear'}
        Interpolation for the registration result (default is 'NearestNeighbor').
    convergence: float
        Threshold for convergence, can make the algorithm very slow (default is
        convergence).
    ignore_affine: bool
        Ignore the affine matrix information extracted from the image header. (default
        is False).
    ignore_header: bool
        Ignore the orientation information and affine matrix information extracted from
        the image header (default is False).

    Notes
    ----------
    The code was mainly taken from nighres.registration.embedded_antsreg and was
    slightly customized. Originally, the module was a port of the CBSTools Java module
    by Pierre-Louis Bazin. The main algorithm is part of the ANTs software by Brian
    Avants and colleagues [1]_. Parameters have been set to values commonly found in
    neuroimaging scripts online, but not necessarily optimal.

    References
    ----------
    .. [1] Avants et al (2008), Symmetric diffeomorphic
       image registration with cross-correlation: evaluating automated labeling
       of elderly and neurodegenerative brain, Med Image Anal. 12(1):26-41
    """
    # check if ants is installed to raise sensible error
    check_installation("antsRegistration")
    check_installation("antsApplyTransforms")

    # make sure that saving related parameters are correct
    _, name_source, ext_source = get_filename(file_source)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    transformed_source_file = os.path.join(
        output_dir, f"{name_source}_ants-def{ext_source}"
    )
    mapping_file = os.path.join(output_dir, f"{name_source}_ants-map{ext_source}")
    inverse_mapping_file = os.path.join(
        output_dir, f"{name_source}_ants-invmap{ext_source}"
    )

    # load and get dimensions and resolution from input images
    source = load_volume(file_source)
    src_affine = source.affine
    nsx = source.header.get_data_shape()[X]
    nsy = source.header.get_data_shape()[Y]
    nsz = source.header.get_data_shape()[Z]
    rsx = source.header.get_zooms()[X]
    rsy = source.header.get_zooms()[Y]
    rsz = source.header.get_zooms()[Z]

    orig_src_aff = source.affine
    orig_src_hdr = source.header

    target = load_volume(file_target)
    trg_affine = target.affine
    ntx = target.header.get_data_shape()[X]
    nty = target.header.get_data_shape()[Y]
    ntz = target.header.get_data_shape()[Z]
    rtx = target.header.get_zooms()[X]
    rty = target.header.get_zooms()[Y]
    rtz = target.header.get_zooms()[Z]

    orig_trg_aff = target.affine
    orig_trg_hdr = target.header

    # in case the affine transformations are not to be trusted: make them equal
    if ignore_affine or ignore_header:
        # create generic affine aligned with the orientation for the source
        new_affine = np.zeros((4, 4))
        if ignore_header:
            new_affine[0][0] = rsx
            new_affine[1][1] = rsy
            new_affine[2][2] = rsz
            new_affine[0][3] = -rsx * nsx / 2.0
            new_affine[1][3] = -rsy * nsy / 2.0
            new_affine[2][3] = -rsz * nsz / 2.0
        else:
            mx = np.argmax(
                np.abs([src_affine[0][0], src_affine[1][0], src_affine[2][0]])
            )
            my = np.argmax(
                np.abs([src_affine[0][1], src_affine[1][1], src_affine[2][1]])
            )
            mz = np.argmax(
                np.abs([src_affine[0][2], src_affine[1][2], src_affine[2][2]])
            )
            new_affine[mx][0] = rsx * np.sign(src_affine[mx][0])
            new_affine[my][1] = rsy * np.sign(src_affine[my][1])
            new_affine[mz][2] = rsz * np.sign(src_affine[mz][2])
            if np.sign(src_affine[mx][0]) < 0:
                new_affine[mx][3] = rsx * nsx / 2.0
            else:
                new_affine[mx][3] = -rsx * nsx / 2.0

            if np.sign(src_affine[my][1]) < 0:
                new_affine[my][3] = rsy * nsy / 2.0
            else:
                new_affine[my][3] = -rsy * nsy / 2.0

            if np.sign(src_affine[mz][2]) < 0:
                new_affine[mz][3] = rsz * nsz / 2.0
            else:
                new_affine[mz][3] = -rsz * nsz / 2.0

        new_affine[3][3] = 1.0

        src_img = nb.Nifti1Image(source.get_fdata(), new_affine, source.header)
        src_img.update_header()
        src_img_file = os.path.join(output_dir, f"tmp_srcimg{ext_source}")
        save_volume(src_img_file, src_img)
        source = load_volume(src_img_file)
        src_affine = source.affine

        # create generic affine aligned with the orientation for the target
        new_affine = np.zeros((4, 4))
        if ignore_header:
            new_affine[0][0] = rtx
            new_affine[1][1] = rty
            new_affine[2][2] = rtz
            new_affine[0][3] = -rtx * ntx / 2.0
            new_affine[1][3] = -rty * nty / 2.0
            new_affine[2][3] = -rtz * ntz / 2.0
        else:
            mx = np.argmax(
                np.abs([trg_affine[0][0], trg_affine[1][0], trg_affine[2][0]])
            )
            my = np.argmax(
                np.abs([trg_affine[0][1], trg_affine[1][1], trg_affine[2][1]])
            )
            mz = np.argmax(
                np.abs([trg_affine[0][2], trg_affine[1][2], trg_affine[2][2]])
            )
            new_affine[mx][0] = rtx * np.sign(trg_affine[mx][0])
            new_affine[my][1] = rty * np.sign(trg_affine[my][1])
            new_affine[mz][2] = rtz * np.sign(trg_affine[mz][2])
            if np.sign(trg_affine[mx][0]) < 0:
                new_affine[mx][3] = rtx * ntx / 2.0
            else:
                new_affine[mx][3] = -rtx * ntx / 2.0

            if np.sign(trg_affine[my][1]) < 0:
                new_affine[my][3] = rty * nty / 2.0
            else:
                new_affine[my][3] = -rty * nty / 2.0

            if np.sign(trg_affine[mz][2]) < 0:
                new_affine[mz][3] = rtz * ntz / 2.0
            else:
                new_affine[mz][3] = -rtz * ntz / 2.0

        new_affine[3][3] = 1.0

        trg_img = nb.Nifti1Image(target.get_fdata(), new_affine, target.header)
        trg_img.update_header()
        trg_img_file = os.path.join(output_dir, f"tmp_trgimg{ext_source}")
        save_volume(trg_img_file, trg_img)
        target = load_volume(trg_img_file)
        trg_affine = target.affine

    # build coordinate mapping matrices and save them to disk
    generate_coordinate_mapping(
        file_source,
        pad=0,
        path_output=output_dir,
        suffix="tmp_srccoord",
        time=False,
        write_output=True,
    )
    generate_coordinate_mapping(
        file_target,
        pad=0,
        path_output=output_dir,
        suffix="tmp_trgcoord",
        time=False,
        write_output=True,
    )
    src_map_file = os.path.join(output_dir, "cmap_tmp_srccoord.nii")
    trg_map_file = os.path.join(output_dir, "cmap_tmp_trgcoord.nii")
    src_map = load_volume(src_map_file)
    trg_map = load_volume(trg_map_file)

    srcfile = source.get_filename()
    trgfile = target.get_filename()

    # run the main ANTS software: here we directly build the command line call
    reg = "antsRegistration --collapse-output-transforms 1 --dimensionality 3 --float 0"
    reg += " --initialize-transforms-per-stage 0 --interpolation Linear"
    reg += " --use-histogram-matching 0"
    reg += " --winsorize-image-intensities [0.001,0.999]"
    reg += f" --initial-moving-transform [{srcfile},{trgfile},1]"
    reg += " --write-composite-transform 0"

    # output with basename antsreg
    output_name = "syn"
    file_out = os.path.join(output_dir, f"{output_name}.nii.gz")
    reg += f" --output [{output_dir}/{output_name},{file_out}]"

    print(f"registering {srcfile} to {trgfile}")

    # figure out the number of scales, going with a factor of two
    scaling_factor = 8
    n_scales = math.ceil(math.log(scaling_factor) / math.log(2.0))
    iter_rigid = str(rigid_iterations)
    iter_affine = str(affine_iterations)
    iter_syn = str(coarse_iterations)
    smooth = str(scaling_factor)
    shrink = str(scaling_factor)
    for n in range(n_scales):
        iter_rigid = f"{iter_rigid}x{rigid_iterations}"
        iter_affine = f"{iter_affine}x{affine_iterations}"
        if n < (n_scales - 1) / 2:
            iter_syn = f"{iter_syn}x{coarse_iterations}"
        elif n < n_scales - 1:
            iter_syn = f"{iter_syn}x{medium_iterations}"
        else:
            iter_syn = f"{iter_syn}x{fine_iterations}"
        smooth = f"{smooth}x{scaling_factor / math.pow(2.0, n + 1)}"
        shrink = f"{shrink}x{math.ceil(scaling_factor / math.pow(2.0, n + 1))}"

    # set parameters for all the different types of transformations
    if run_rigid is True:
        reg += " --transform Rigid[0.1]"
        if cost_function == "CrossCorrelation":
            reg += f" --metric CC[{trgfile},{srcfile},1.000,5,Random,0.3]"
        else:
            reg += f" --metric MI[{trgfile},{srcfile},1.000,32,Random,0.3]"

        reg += f" --convergence [{iter_rigid},{convergence},10]"
        reg += f" --smoothing-sigmas {smooth}"
        reg += f" --shrink-factors {shrink}"

    if run_affine is True:
        reg += " --transform Affine[0.1]"
        if cost_function == "CrossCorrelation":
            reg += f" --metric CC[{trgfile},{srcfile},1.000,5,Random,0.3]"
        else:
            reg += f" --metric MI[{trgfile},{srcfile},1.000,32,Random,0.3]"

        reg += f" --convergence [{iter_affine},{convergence},10]"
        reg += f" --smoothing-sigmas {smooth}"
        reg += f" --shrink-factors {shrink}"

    if run_syn is True:
        # Regularization preset for the SyN deformation
        # syn_param = [0.2, 1.0, 0.0] # Regularization Low
        syn_param = [0.2, 3.0, 0.0]  # Regularization Medium
        # syn_param = [0.2, 4.0, 3.0] # Regularization High

        reg += f" --transform SyN{syn_param}"
        if cost_function == "CrossCorrelation":
            reg += f" --metric CC[{trgfile},{srcfile},1.000,5,Random,0.3]"
        else:
            reg += f" --metric MI[{trgfile},{srcfile},1.000,32,Random,0.3]"

        reg += f" --convergence [{iter_syn},{convergence},5]"
        reg += f" --smoothing-sigmas {smooth}"
        reg += f" --shrink-factors {shrink}"

    if run_rigid is False and run_affine is False and run_syn is False:
        reg += " --transform Rigid[0.1]"
        reg += f" --metric CC[{trgfile},{srcfile},1.000,5,Random,0.3]"
        reg += " --convergence [0,1.0,2]"
        reg += " --smoothing-sigmas 0.0"
        reg += " --shrink-factors 1"

    # run the ANTs command directly
    execute_command(reg)

    # output file names
    results = sorted(glob(os.path.join(output_dir, f"{output_name}*")))
    print(results)
    forward = []
    flag = []
    for res in results:
        if res.endswith("GenericAffine.mat"):
            forward.append(res)
            flag.append(False)
        elif res.endswith("Warp.nii.gz") and not res.endswith("InverseWarp.nii.gz"):
            forward.append(res)
            flag.append(False)

    # convert transforms to coordinate mappings in cbs format
    inverse = []
    linear = []
    for res in results[::-1]:
        if res.endswith("GenericAffine.mat"):
            inverse.append(res)
            linear.append(True)
        elif res.endswith("InverseWarp.nii.gz"):
            inverse.append(res)
            linear.append(False)

    at = "antsApplyTransforms --dimensionality 3 --input-image-type 0"
    at += f" --input {source.get_filename()}"
    at += f" --reference-image {target.get_filename()}"
    at += f" --interpolation {interpolation}"
    for idx, transform in enumerate(forward):
        if flag[idx]:
            at += f" --transform [{transform}, 1]"
        else:
            at += f" --transform [{transform}, 0]"
        at += f" --output {transformed_source_file}"

    # run
    execute_command(at)

    # Create coordinate mappings
    src_at = "antsApplyTransforms --dimensionality 3 --input-image-type 3"
    src_at += f" --input {src_map.get_filename()}"
    src_at += f" --reference-image {target.get_filename()}"
    src_at += " --interpolation Linear"
    for idx, transform in enumerate(forward):
        if flag[idx]:
            src_at += f" --transform [{transform}, 1]"
        else:
            src_at += f" --transform [{transform}, 0]"
    src_at += f" --output {mapping_file}"

    # run
    execute_command(src_at)

    trg_at = "antsApplyTransforms --dimensionality 3 --input-image-type 3"
    trg_at += f" --input {trg_map.get_filename()}"
    trg_at += f" --reference-image {source.get_filename()}"
    trg_at += " --interpolation Linear"
    for idx, transform in enumerate(inverse):
        if linear[idx]:
            trg_at += f" --transform [{transform}, 1]"
        else:
            trg_at += f" --transform [{transform}, 0]"
    trg_at += f" --output {inverse_mapping_file}"

    # run
    execute_command(trg_at)

    # clean-up intermediate files
    if os.path.exists(file_out):
        os.remove(file_out)
    if os.path.exists(src_map_file):
        os.remove(src_map_file)
    if os.path.exists(trg_map_file):
        os.remove(trg_map_file)
    if ignore_affine or ignore_header:
        if os.path.exists(src_img_file):
            os.remove(src_img_file)
        if os.path.exists(trg_img_file):
            os.remove(trg_img_file)

    for name in forward:
        if os.path.exists(name):
            os.remove(name)
    for name in inverse:
        if os.path.exists(name):
            os.remove(name)

    # if ignoring header and/or affine, must paste back the correct headers
    if ignore_affine or ignore_header:
        mapping = load_volume(mapping_file)
        save_volume(
            mapping_file,
            nb.Nifti1Image(mapping.get_fdata(), orig_trg_aff, orig_trg_hdr),
        )
        inverse = load_volume(inverse_mapping_file)
        save_volume(
            inverse_mapping_file,
            nb.Nifti1Image(inverse.get_fdata(), orig_src_aff, orig_src_hdr),
        )
        trans = load_volume(transformed_source_file)
        save_volume(
            transformed_source_file,
            nb.Nifti1Image(trans.get_fdata(), orig_trg_aff, orig_trg_hdr),
        )
    print("done")


def recursive_bbr(
    file_ref,
    lh_white,
    rh_white,
    lh_pial,
    rh_pial,
    subjects_dir,
    file_mask=None,
    min_vox=10,
    min_vtx=500,
    cuboid=True,
    tetrahedra=False,
    nn_smooth=0.5,
    registration_mode="syty",
    reverse_contrast=True,
    cost_method="GreveFischl",
    contrast_method="gradient",
):
    """This function performs recurive BBR to register a functional data set with
    surfaces generated from a separate whole-brain anatomy.

    parameter mode (explanation from one of Tim vanMourik's scripts):
    `r' for rotation around all axes, `rx', `ry', or `rz' for a single rotation around
    the respective axis and for example `rxrz' for a rotation around the x-axis and the
    z-axis. In the same way the scaling and translation can be modified. For example the
    combination `rystytz' would optimise for a rotation around the y-axis, a scaling in
    all direction and a translation along the y-axis and along the z-axis. The default
    is 'rst' for all transformations. The script needs an installation of freesurfer.

    Dependencies (should be added to the matlab search path):
    OpenFmriAnalysis: https://github.com/TimVanMourik/OpenFmriAnalysis

    Parameters
    ----------
    file_ref : str
        File name of reference volume.
    lh_white : str
        File name of left white surface.
    rh_white : str
        File name of right white surface.
    lh_pial : str
        File name of left pial surface.
    rh_pial : str
        File name of right pial surface.
    subjects_dir : str
        Subjects directory.
    file_mask : str, optional
        File name of mask volume, by default None.
    min_vox : int, optional
        Minimum number of voxels, by default 10.
    min_vtx : int, optional
        Minimum number of vertices, by default 500.
    cuboid : bool, optional
        Divide into smaller cubes, by default True.
    tetrahedra : bool, optional
        Divide into smaller tetrahedra, by default False.
    nn_smooth : float, optional
        Smoothing parameter, by default 0.5.
    registration_mode : str, optional
        Registration mode as explained further above, by default "syty".
    reverse_contrast : bool, optional
        If inner surface is darker than outer surface, by default True.
    cost_method : str, optional
        Can be GreveFischl, sum, sumOfSquares, by default "GreveFischl".
    contrast_method : str, optional
        Can be gradient, average, fixedDistance, relativeDistance, extrema, by default
        "gradient".
    """
    # make output folder
    if not exists(subjects_dir):
        os.makedirs(subjects_dir)

    path_input = join(subjects_dir, "input")
    if not exists(path_input):
        os.makedirs(path_input)

    # copy input files to input folder
    sh.copyfile(lh_white, join(path_input, "lh.white"))
    sh.copyfile(rh_white, join(path_input, "rh.white"))
    sh.copyfile(lh_pial, join(path_input, "lh.pial"))
    sh.copyfile(rh_pial, join(path_input, "rh.pial"))
    sh.copyfile(file_ref, join(path_input, "reference_volume.nii"))
    if file_mask is not None:
        sh.copyfile(file_mask, join(path_input, "mask.nii"))

    # filenames
    input_surf_ras = "boundaries_in_ras.mat"
    input_surf_vox = "boundaries_in_vox.mat"
    output_surf_ras = "boundaries_out_ras.mat"
    output_surf_vox = "boundaries_out_vox.mat"
    output_cmap = "cmap_rBBR.nii"

    lh_white = join(path_input, "lh.white")
    rh_white = join(path_input, "rh.white")
    lh_pial = join(path_input, "lh.pial")
    rh_pial = join(path_input, "rh.pial")
    ref_vol = join(path_input, "reference_volume.nii")
    in_surf_mat_ras = join(path_input, input_surf_ras)
    in_surf_mat_vox = join(path_input, input_surf_vox)
    out_surf_mat_ras = join(subjects_dir, output_surf_ras)
    out_surf_mat_vox = join(subjects_dir, output_surf_vox)
    input_cfg = join(path_input, "cfg.mat")

    # get affine vox2ras-tkr and ras2vox-tkr transformation to reference volume
    vox2ras_tkr, ras2vox_tkr = read_vox2ras_tkr(ref_vol)

    # surf2mat
    matlab = MatlabCommand(
        "ft_surf2mat", lh_white, rh_white, lh_pial, rh_pial, in_surf_mat_ras
    )
    matlab.run()

    # ras2vox
    data = loadmat(in_surf_mat_ras)

    # apply ras2vox transformation to vertices
    for i in range(2):
        data["wSurface"][0][i] = apply_affine_chunked(
            ras2vox_tkr, data["wSurface"][0][i]
        )
        data["pSurface"][0][i] = apply_affine_chunked(
            ras2vox_tkr, data["pSurface"][0][i]
        )

    # save surfaces in voxel space
    savemat(in_surf_mat_vox, data)

    # get configuration mat-file for rBBR
    cfg = {}
    cfg["subjects_dir"] = subjects_dir
    cfg["reference_volume"] = join("input", "reference_volume.nii")
    cfg["input_surf"] = join("input", input_surf_vox)
    cfg["mask"] = join("input", "mask.nii") if len(file_mask) else ""
    cfg["min_vox"] = float(min_vox)
    cfg["min_vtx"] = float(min_vtx)
    cfg["cuboid"] = cuboid
    cfg["tetrahedra"] = tetrahedra
    cfg["nn_smooth"] = float(nn_smooth)
    cfg["registration_mode"] = registration_mode
    cfg["reverse_contrast"] = reverse_contrast
    cfg["cost_method"] = cost_method
    cfg["contrast_method"] = contrast_method
    cfg["output_surf"] = output_surf_vox
    cfg["output_cmap"] = output_cmap

    savemat(input_cfg, cfg)

    # run rBBR
    matlab = MatlabCommand("ft_rBBR", input_cfg)
    matlab.run()

    # vox2ras
    data = loadmat(out_surf_mat_vox)

    # apply ras2vox transformation to output
    for i in range(2):
        data["wSurface"][0][i] = data["wSurface"][0][i][:, 0:3]
        data["pSurface"][0][i] = data["pSurface"][0][i][:, 0:3]

        data["wSurface"][0][i] = apply_affine_chunked(
            vox2ras_tkr, data["wSurface"][0][i]
        )
        data["pSurface"][0][i] = apply_affine_chunked(
            vox2ras_tkr, data["pSurface"][0][i]
        )

    # save matfile
    savemat(out_surf_mat_ras, data)

    # save surfaces in freesurfer format
    write_geometry(
        join(subjects_dir, "lh.white"), data["wSurface"][0][0], data["faceData"][0][0]
    )
    write_geometry(
        join(subjects_dir, "rh.white"), data["wSurface"][0][1], data["faceData"][0][1]
    )
    write_geometry(
        join(subjects_dir, "lh.pial"), data["pSurface"][0][0], data["faceData"][0][0]
    )
    write_geometry(
        join(subjects_dir, "rh.pial"), data["pSurface"][0][1], data["faceData"][0][1]
    )
