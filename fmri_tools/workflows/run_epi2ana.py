# -*- coding: utf-8 -*-
"""EPI <-> ANA registration

The purpose of the following script is to compute the deformation field for the
registration between antomy and EPI in native space. The mask should be in ana
space. Optionally, a second mask can be applied which must be in orig space. The
script consists of the following steps:
    1. enhance brainmask if second mask is given (brain.finalsurf.mgz)
    2. set output folder structure
    3. scanner transform t1 <-> epi
    4. n4 correction epi
    5. clean ana (remove ceiling and normalise)
    6. mask t1 and epi
    7. antsreg
    8. clean deformations
    9. expand deformations
    10. apply deformation

The script needs an installation of freesurfer and ants.

"""

import os
import shutil as sh
from argparse import SUPPRESS, ArgumentParser

import nibabel as nb

from ..registration.cmap import clean_coordinate_mapping, expand_coordinate_mapping
from ..registration.nonrigid import embedded_antsreg
from ..registration.transform import apply_coordinate_mapping
from ..segmentation.mask import clean_ana, mask_ana, mask_epi
from ..segmentation.skullstrip import skullstrip_refined
from ..segmentation.vol import remove_bias_ants

# cmap parameters
CLEAN_CMAP = True
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

# compress output files
COMPRESS_FILE = False


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
        "--epi",
        dest="file_epi",
        type=str,
        help="File name of epi volume.",
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
        "--mask2",
        dest="file_mask2",
        type=int,
        help=(
            "Second image used as brain mask (brain.finalsurf) in freesurfer space (default: %(default)s)."
        ),
        default=None,
    )
    optional.add_argument(
        "--clean",
        dest="cleanup",
        type=int,
        help=("Remove intermediate files (default: %(default)s)."),
        default=False,
    )

    return parser


def epi2ana_workflow(file_epi, file_t1, file_mask, path_output, file_mask2, cleanup):
    """Workflow for epi2ana registration.

    Parameters
    ----------
    file_epi : str
        File name of mean epi volume.
    file_t1 : str
        File name of MP2RAGE T1 image.
    file_mask : str
        File name of MP2RAGE skullstrip mask.
    path_output : str
        Directory where output is saved.
    file_mask2 : str
        Optional mask (brain.finalsurf) in freesurfer space.
    cleanup : bool
        Remove intermediate files.
    """
    # enhance brainmask
    if file_mask2 is not None:
        file_mask = skullstrip_refined(file_mask, file_mask2)

    # set folder structure
    path_temp = os.path.join(path_output, "temp")
    path_epi = os.path.join(path_temp, "epi")
    path_t1 = os.path.join(path_temp, "t1")
    path_syn = os.path.join(path_temp, "syn")

    if not os.path.exists(path_output):
        os.makedirs(path_output)

    if not os.path.exists(path_temp):
        os.makedirs(path_temp)

    if not os.path.exists(path_epi):
        os.makedirs(path_epi)

    if not os.path.exists(path_t1):
        os.makedirs(path_t1)

    if not os.path.exists(path_syn):
        os.makedirs(path_syn)

    # copy input files
    sh.copyfile(file_epi, os.path.join(path_epi, "epi.nii"))
    sh.copyfile(file_t1, os.path.join(path_t1, "T1.nii"))
    sh.copyfile(file_mask, os.path.join(path_t1, "mask.nii"))

    # bias field correction to epi
    remove_bias_ants(
        os.path.join(path_epi, "epi.nii"), os.path.join(path_epi, "bepi.nii")
    )

    # clean ana
    clean_ana(os.path.join(path_t1, "T1.nii"), 1000.0, 4095.0, overwrite=True)

    # mask t1 and epi
    mask_ana(
        os.path.join(path_t1, "T1.nii"),
        os.path.join(path_t1, "mask.nii"),
        background_bright=False,
    )
    mask_epi(
        os.path.join(path_epi, "bepi.nii"),
        os.path.join(path_t1, "pT1.nii"),
        os.path.join(path_t1, "mask.nii"),
        NITER_MASK,
        SIGMA_MASK,
    )

    # syn
    embedded_antsreg(
        os.path.join(path_t1, "pT1.nii"),  # source image
        os.path.join(path_epi, "pbepi.nii"),  # target image
        path_syn,
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

    # output file names
    fout_epi2ana = os.path.join(path_output, "epi2ana.nii")
    fout_epi2ana_example = os.path.join(path_output, "epi2ana_example.nii")
    fout_ana2epi = os.path.join(path_output, "ana2epi.nii")
    fout_ana2epi_example = os.path.join(path_output, "ana2epi_example.nii")
    fout_mask = os.path.join(path_output, "epi2ana_mask.nii")

    # rename final deformations
    os.rename(os.path.join(path_syn, "pT1_ants-map.nii"), fout_ana2epi)
    os.rename(os.path.join(path_syn, "pT1_ants-invmap.nii"), fout_epi2ana)

    # clean deformation
    if CLEAN_CMAP:
        epi2ana_cleaned = clean_coordinate_mapping(
            fout_ana2epi,
            fout_epi2ana,
            overwrite_file=True,
            save_mask=False,
        )

        # write mask
        nb.save(epi2ana_cleaned["mask"], fout_mask)

    # expand deformation
    if EXPAND_CMAP:
        _ = expand_coordinate_mapping(
            fout_ana2epi,
            path_output,
            name_output="ana2epi",
            write_output=True,
        )

        _ = expand_coordinate_mapping(
            fout_epi2ana,
            path_output,
            name_output="epi2ana",
            write_output=True,
        )

    # apply deformation
    # ana -> epi
    apply_coordinate_mapping(
        file_t1,  # input
        fout_ana2epi,  # cmap
        fout_ana2epi_example,
        interpolation="linear",  # nearest or linear
    )

    # epi -> ana
    apply_coordinate_mapping(
        file_epi,  # input
        fout_epi2ana,  # cmap
        fout_epi2ana_example,
        interpolation="linear",  # nearest or linear
    )

    if COMPRESS_FILE:
        nb.save(nb.load(fout_ana2epi), f"{fout_ana2epi}.gz")
        nb.save(nb.load(fout_epi2ana), f"{fout_epi2ana}.gz")
        nb.save(nb.load(fout_ana2epi_example), f"{fout_ana2epi_example}.gz")
        nb.save(nb.load(fout_epi2ana_example), f"{fout_epi2ana_example}.gz")
        nb.save(nb.load(fout_mask), f"{fout_mask}.gz")

        os.remove(fout_ana2epi)
        os.remove(fout_epi2ana)
        os.remove(fout_ana2epi_example)
        os.remove(fout_epi2ana_example)
        os.remove(fout_mask)

    # clean intermediate files
    if cleanup:
        sh.rmtree(path_temp, ignore_errors=True)


def _main(argv=None):
    """EPI2ANA registration workflow."""
    options = _get_parser().parse_args(argv)
    kwargs = vars(options)
    epi2ana_workflow(**kwargs)


if __name__ == "__main__":
    _main()
