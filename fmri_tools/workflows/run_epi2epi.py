# -*- coding: utf-8 -*-
"""EPI <-> EPI

The purpose of the following script is to compute the deformation field for the
registration between different epi time series. The script consists of the
following steps:
    1. set output folder structure
    2. n4 correction epi
    3. clean ana (remove ceiling and normalise)
    4. mask epi
    5. antsreg or flirt
    6. expand coordinate mapping
    7. apply deformations

The script needs an installation of freesurfer and ants (or fsl).

"""

import os
import shutil as sh
import subprocess
from argparse import SUPPRESS, ArgumentParser

from ..registration.cmap import expand_coordinate_mapping, generate_coordinate_mapping
from ..registration.fsl import apply_flirt, flirt
from ..registration.nonrigid import embedded_antsreg
from ..registration.transform import apply_coordinate_mapping
from ..segmentation.mask import clean_ana, mask_ana, mask_epi
from ..segmentation.vol import remove_bias_ants

# cmap parameters
EXPAND_CMAP = True

# parameters for epi skullstrip
NITER_MASK = 3
SIGMA_MASK = 3

# parameters for syn
RUN_RIGID = True
RIGID_ITERATIONS = 1000
RUN_AFFINE = False
AFFINE_ITERATIONS = 1000
RUN_SYN = True
COARSE_ITERATIONS = 50
MEDIUM_ITERATIONS = 150
FINE_ITERATIONS = 100
COST_FUNCTION = "CrossCorrelation"
INTERPOLATION = "Linear"


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
        "--source",
        dest="file_epi_source",
        type=str,
        help="File name of epi source image.",
        required=True,
    )
    required.add_argument(
        "--target",
        dest="file_epi_target",
        type=str,
        help="File name of epi target image.",
        required=True,
    )
    required.add_argument(
        "--t1",
        dest="file_t1",
        type=str,
        help="File name of MP2RAGE image.",
        required=True,
    )
    required.add_argument(
        "--mask",
        dest="file_mask",
        type=str,
        help="File name of MP2RAGE skullstrip mask.",
        required=True,
    )
    required.add_argument(
        "--out",
        dest="path_output",
        type=str,
        help="Directory where output is saved.",
        required=True,
    )
    optional.add_argument(
        "--syn",
        dest="syn",
        type=int,
        help="Apply Syn if True else Flirt (default: %(default)s).",
        default=True,
    )
    optional.add_argument(
        "--clean",
        dest="cleanup",
        type=int,
        help=("Remove intermediate files (default: %(default)s)."),
        default=False,
    )

    return parser


def epi2epi_workflow(
    file_epi_source,
    file_epi_target,
    file_t1,
    file_mask,
    path_output,
    syn,
    cleanup,
):
    """Workflow for epi2epi registration.

    Parameters
    ----------
    file_epi_source : str
        File name of epi source image.
    file_epi_target : str
        File name of epi target image.
    file_t1 : str
        File name of MP2RAGE T1 image.
    file_mask : str
        File name of MP2RAGE skullstrip mask.
    path_output : str
        Directory where output is saved.
    syn:: bool
        Apply syn (ants) if true else flirt (fsl).
    cleanup : bool
        Remove intermediate files.
    """
    # set folder structure
    path_temp = os.path.join(path_output, "temp")
    path_epi_source = os.path.join(path_temp, "epi_source")
    path_epi_target = os.path.join(path_temp, "epi_target")
    path_t1_source = os.path.join(path_temp, "t1_source")
    path_t1_target = os.path.join(path_temp, "t1_target")
    path_reg = os.path.join(path_temp, "reg")

    if not os.path.exists(path_output):
        os.makedirs(path_output)

    if not os.path.exists(path_temp):
        os.makedirs(path_temp)

    if not os.path.exists(path_epi_source):
        os.makedirs(path_epi_source)

    if not os.path.exists(path_epi_target):
        os.makedirs(path_epi_target)

    if not os.path.exists(path_t1_source):
        os.makedirs(path_t1_source)

    if not os.path.exists(path_t1_target):
        os.makedirs(path_t1_target)

    if not os.path.exists(path_reg):
        os.makedirs(path_reg)

    path_t1 = [path_t1_source, path_t1_target]
    path_epi = [path_epi_source, path_epi_target]

    # copy input files
    sh.copyfile(file_epi_source, os.path.join(path_epi_source, "epi.nii"))
    sh.copyfile(file_epi_target, os.path.join(path_epi_target, "epi.nii"))
    sh.copyfile(file_t1, os.path.join(path_t1_source, "T1.nii"))
    sh.copyfile(file_mask, os.path.join(path_t1_source, "mask.nii"))
    sh.copyfile(file_t1, os.path.join(path_t1_target, "T1.nii"))
    sh.copyfile(file_mask, os.path.join(path_t1_target, "mask.nii"))

    # bias field correction to epi
    for i, _ in enumerate(path_epi):
        remove_bias_ants(
            os.path.join(path_epi[i], "epi.nii"), os.path.join(path_epi[i], "bepi.nii")
        )

    # clean ana
    for i, _ in enumerate(path_t1):
        clean_ana(os.path.join(path_t1[i], "T1.nii"), 1000.0, 4095.0, overwrite=True)

    # mask t1 and epi
    for i, _ in enumerate(path_t1):
        mask_ana(
            os.path.join(path_t1[i], "T1.nii"),
            os.path.join(path_t1[i], "mask.nii"),
            background_bright=False,
        )

    for i, _ in enumerate(path_epi):
        mask_epi(
            os.path.join(path_epi[i], "bepi.nii"),
            os.path.join(path_t1[i], "pT1.nii"),
            os.path.join(path_t1[i], "mask.nii"),
            NITER_MASK,
            SIGMA_MASK,
        )

    if syn:
        print("Apply SyN registration using ANTS.")
        _apply_syn(
            os.path.join(path_epi_target, "pbepi.nii"),  # source image
            os.path.join(path_epi_source, "pbepi.nii"),  # target image
            path_output,
        )
    else:
        print("Apply flirt registration using FSL.")
        _apply_flirt(
            os.path.join(path_epi_target, "pbepi.nii"),  # source image
            os.path.join(path_epi_source, "pbepi.nii"),  # target image
            path_output,
        )

    # expand deformation
    if EXPAND_CMAP:
        _ = expand_coordinate_mapping(
            os.path.join(path_output, "source2target.nii.gz"),
            path_output,
            name_output="source2target",
            write_output=True,
        )

        _ = expand_coordinate_mapping(
            os.path.join(path_output, "target2source.nii.gz"),
            path_output,
            name_output="target2source",
            write_output=True,
        )

    # apply deformation
    # source -> target
    apply_coordinate_mapping(
        file_epi_source,  # input
        os.path.join(path_output, "source2target.nii.gz"),  # cmap
        os.path.join(path_output, "source2target_example.nii.gz"),
        interpolation="linear",  # nearest or linear
    )

    # target -> source
    apply_coordinate_mapping(
        file_epi_target,  # input
        os.path.join(path_output, "target2source.nii.gz"),  # cmap
        os.path.join(path_output, "target2source_example.nii.gz"),
        interpolation="linear",  # nearest or linear
    )

    # clean intermediate files
    if cleanup:
        sh.rmtree(path_temp, ignore_errors=True)


def _apply_syn(source_in, target_in, dir_out):
    """Apply SyN registration.

    Parameters
    ----------
    source_in : str
        File name of source image.
    target_in : str
        File name of target image.
    dir_out : str
        Directory where output is saved.
    """
    embedded_antsreg(
        source_in,
        target_in,
        dir_out,
        RUN_RIGID,  # whether or not to run a rigid registration first
        RIGID_ITERATIONS,  # number of iterations in the rigid step
        RUN_AFFINE,  # whether or not to run an affine registration first
        AFFINE_ITERATIONS,  # number of iterations in the affine step
        RUN_SYN,  # whether or not to run a SyN registration
        COARSE_ITERATIONS,  # number of iterations at the coarse level
        MEDIUM_ITERATIONS,  # number of iterations at the medium level
        FINE_ITERATIONS,  # number of iterations at the fine level
        COST_FUNCTION,  # CrossCorrelation or MutualInformation
        INTERPOLATION,  # interpolation for registration result (NeareastNeighbor or Linear)
        convergence=1e-6,  # threshold for convergence (can make algorithm very slow)
        ignore_affine=False,  # ignore the affine matrix information extracted from the image header
        ignore_header=False,  # ignore the orientation information and affine matrix information extracted from the image header
    )

    # rename final deformations
    os.rename(
        os.path.join(dir_out, "pbepi_ants-map.nii.gz"),
        os.path.join(dir_out, "target2source.nii.gz"),
    )
    os.rename(
        os.path.join(dir_out, "pbepi_ants-invmap.nii.gz"),
        os.path.join(dir_out, "source2target.nii.gz"),
    )


def _apply_flirt(source_in, target_in, dir_out):
    """Apply Flirt registration.

    Parameters
    ----------
    source_in : str
        File name of source image.
    target_in : str
        File name of target image.
    dir_out : str
        Directory where output is saved.
    """
    os.chdir(dir_out)
    flirt(
        source_in,
        target_in,
        os.path.join(dir_out, "flirt.nii"),
        os.path.join(dir_out, "flirt_matrix.mat"),
        cost_func="corratio",
        interp_method="trilinear",
    )

    # invert matrix
    command = "convert_xfm"
    command += f" -omat {os.path.join(dir_out, 'flirt_inv_matrix.mat')}"
    command += f" -inverse {os.path.join(dir_out, 'flirt_matrix.mat')}"

    print("Execute: " + command)
    try:
        subprocess.run([command], shell=True, check=False)
    except subprocess.CalledProcessError:
        print("Execuation failed!")

    # get cmap
    generate_coordinate_mapping(
        target_in,
        pad=0,
        path_output=dir_out,
        suffix="target",
        time=False,
        write_output=True,
    )

    generate_coordinate_mapping(
        source_in,
        pad=0,
        path_output=dir_out,
        suffix="source",
        time=False,
        write_output=True,
    )

    # apply flirt to cmap
    apply_flirt(
        os.path.join(dir_out, "cmap_target.nii"),
        source_in,
        os.path.join(dir_out, "flirt_matrix.mat"),
        os.path.join(dir_out, "target2source.nii.gz"),
        "trilinear",
    )

    apply_flirt(
        os.path.join(dir_out, "cmap_source.nii"),
        target_in,
        os.path.join(dir_out, "flirt_inv_matrix.mat"),
        os.path.join(dir_out, "source2target.nii.gz"),
        "trilinear",
    )


def _main(argv=None):
    """EPI2EPI registration workflow."""
    options = _get_parser().parse_args(argv)
    kwargs = vars(options)
    epi2epi_workflow(**kwargs)


if __name__ == "__main__":
    _main()
