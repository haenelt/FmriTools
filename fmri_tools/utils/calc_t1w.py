# -*- coding: utf-8 -*-

# python standard library inputs
import argparse

# external inputs
import numpy as np
import nibabel as nb
from nipype.interfaces.ants import N4BiasFieldCorrection


def calc_t1w(file_bold, file_vaso, file_out, apply_bias):
    """Compute a T1-weighted image from nulled and not-nulled time series. Be
    aware that the computation takes all time points. Therefore, volumes not in
    steady-state should be removed before running this function.

    Parameters
    ----------
    file_bold : str
        File name of not-nulled time series.
    file_vaso : str
        File name of nulled time series.
    file_out : str
        File name of T1-weighted image.
    apply_bias : bool
        Apply a bias field correction.

    Returns
    -------
    None.

    """

    print("Compute T1-weighted image from")
    print("BOLD: "+file_bold)
    print("VASO: "+file_vaso)

    # load data
    bold = nb.load(file_bold)
    vaso = nb.load(file_vaso)
    arr_bold = bold.get_fdata()
    arr_vaso = vaso.get_fdata()

    bold_dim = bold.header["dim"]
    nx = bold_dim[1]
    ny = bold_dim[2]
    nz = bold_dim[3]
    nt = bold_dim[4]

    arr_combined = np.zeros((nx, ny, nz, 2*nt))
    arr_combined[:, :, :, :nt] = arr_bold
    arr_combined[:, :, :, nt:] = arr_vaso

    arr_abs = np.abs(np.mean(arr_combined, axis=3))
    arr_std = np.std(arr_combined, axis=3)
    arr_std[arr_std == 0] = 1
    arr_std[~np.isfinite(arr_std)] = 1
    arr_res = arr_abs / arr_std

    # write output
    bold.header["dim"][0] = 3
    bold.header["dim"][4] = 1
    output = nb.Nifti1Image(arr_res, bold.affine, bold.header)
    nb.save(output, file_out)

    # bias field correction to epi
    if apply_bias:
        print("Apply bias field correction")
        n4 = N4BiasFieldCorrection()
        n4.inputs.dimension = 3
        n4.inputs.input_image = file_out
        n4.inputs.output_image = file_out
        n4.run()


if __name__ == "__main__":

    # description
    parser_description = "This program computes a T1-weighted image from a" \
                         "nulled and not-nulled time series. All time points" \
                         "are used, i.e., non steady-state volumes should be" \
                         "removed from the time series before calling this " \
                         "program."

    vaso_help = "input file name of nulled time series."
    bold_help = "input file name of not-nulled time series."
    out_help = "output file name."
    bias_help = "apply no bias field correction (optional)."

    # parse arguments from command line
    parser = argparse.ArgumentParser(description=parser_description)
    parser.add_argument('-v', '--vaso', type=str, help=vaso_help)
    parser.add_argument('-b', '--bold', type=str, help=bold_help)
    parser.add_argument('-o', '--out', type=str, help=out_help)
    parser.add_argument('-n', '--no_bias', action='store_false',
                        help=bias_help)
    args = parser.parse_args()

    # run function
    calc_t1w(args.bold, args.vaso, args.out, args.no_bias)

