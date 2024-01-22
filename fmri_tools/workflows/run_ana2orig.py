# -*- coding: utf-8 -*-
"""ANA <-> ORIG registration

The purpose of the following script is to compute the deformation field for the
registration between antomy in conformed freesurfer space and native anatomical
space. The transformation can be computed from the headers of source and target
images. The script needs an installation of freesurfer.

"""

import os
import shutil as sh
from argparse import SUPPRESS, ArgumentParser

from ..io.filename import get_filename
from ..io.vol import mri_convert
from ..registration.transform import apply_coordinate_mapping, scanner_transform


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
        "--t1",
        dest="file_t1",
        type=str,
        help="File name of MP2RAGE image.",
        required=True,
    )
    required.add_argument(
        "--orig",
        dest="file_orig",
        type=str,
        help="File name of freesurfer imag (orig.mgz).",
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
        "--clean",
        dest="cleanup",
        type=int,
        help=("Remove intermediate files (default: %(default)s)."),
        default=False,
    )

    return parser


def ana2orig_workflow(file_t1, file_orig, path_output, cleanup):
    """Workflow for ana2orig transformation.

    Parameters
    ----------
    file_t1 : str
        File name of MP2RAGE image.
    file_orig : str
        File name of corresponding freesurfer image (orig.mgz).
    path_output : str
        Directory where output is saved.
    cleanup : bool
        Remove intermediate files.
    """
    # set folder structure
    path_temp = os.path.join(path_output, "temp")

    if not os.path.exists(path_output):
        os.makedirs(path_output)

    if not os.path.exists(path_temp):
        os.makedirs(path_temp)

    # get filenames
    _, _, ext_orig = get_filename(file_orig)
    _, _, ext_t1 = get_filename(file_t1)

    # copy input files
    sh.copyfile(file_orig, os.path.join(path_temp, "orig" + ext_orig))
    sh.copyfile(file_t1, os.path.join(path_temp, "T1" + ext_t1))

    # convert to nifti
    if ext_orig != ".nii":
        mri_convert(
            os.path.join(path_temp, "orig" + ext_orig),
            os.path.join(path_temp, "orig.nii"),
        )

    if ext_t1 != ".nii":
        mri_convert(
            os.path.join(path_temp, "T1" + ext_t1), os.path.join(path_temp, "T1.nii")
        )

    # scanner transformation
    scanner_transform(
        os.path.join(path_temp, "orig.nii"),
        os.path.join(path_temp, "T1.nii"),
        path_temp,
        True,
    )
    scanner_transform(
        os.path.join(path_temp, "T1.nii"),
        os.path.join(path_temp, "orig.nii"),
        path_temp,
        True,
    )

    # get output
    os.rename(
        os.path.join(path_temp, "orig_2_T1_scanner.nii.gz"),
        os.path.join(path_output, "orig2T1.nii.gz"),
    )
    os.rename(
        os.path.join(path_temp, "T1_2_orig_scanner.nii.gz"),
        os.path.join(path_output, "T12orig.nii.gz"),
    )

    # apply deformation
    # ana -> epi
    apply_coordinate_mapping(
        os.path.join(path_temp, "orig.nii"),  # input
        os.path.join(path_output, "orig2T1.nii.gz"),  # cmap
        os.path.join(path_output, "orig2T1_example.nii.gz"),
        interpolation="linear",  # nearest or linear
    )

    # epi -> ana
    apply_coordinate_mapping(
        os.path.join(path_temp, "T1.nii"),  # input
        os.path.join(path_output, "T12orig.nii.gz"),  # cmap
        os.path.join(path_output, "T12orig_example.nii.gz"),
        interpolation="linear",  # nearest or linear
    )

    # clean intermediate files
    if cleanup:
        sh.rmtree(path_temp, ignore_errors=True)


def _main(argv=None):
    """ANA2ORIG transformation workflow."""
    options = _get_parser().parse_args(argv)
    kwargs = vars(options)
    ana2orig_workflow(**kwargs)


if __name__ == "__main__":
    _main()
