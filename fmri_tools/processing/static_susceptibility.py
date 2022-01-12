# -*- coding: utf-8 -*-

# python standard library inputs
import os
import argparse

# external inputs
import numpy as np
import nibabel as nb
from nipype.interfaces.ants import N4BiasFieldCorrection

# local inputs
from fmri_tools.io.get_filename import get_filename
from fmri_tools.matlab import MatlabCommand


def static_susceptibility(file_in, tr, cutoff_highpass=270, mode="mean",
                          apply_bias=True):
    """Computation of the temporal mean of an EPI time series. A temporal
    baseline correction and a bias field correction on the resulting temporal
    mean is performed. The resulting image can be used to show locations with
    pronounced static susceptibility effects (veins) in functional time series.

    Parameters
    ----------
    file_in : str
        File name of input time series.
    tr : float
        Volume repetition time (TR) in s.
    cutoff_highpass : float, optional
        Cutoff time in s for baseline correction. Not applied if cutoff is set
        to zero.
    mode : str, optional
        Mode for temporal average (mean or median).
    apply_bias : bool, optional
        Apply bias field correction on temporal mean.

    Returns
    -------
    None.

    """

    if mode not in ["mean", "median"]:
        raise ValueError("Average mode must be either mean or median!")

    # read time series
    data = nb.load(file_in)
    arr = data.get_fdata()

    # output folder structure
    path_file, name_file, ext_file = get_filename(file_in)
    path_output = os.path.join(path_file, "ssw", "native")
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # baseline correction
    if cutoff_highpass != 0:
        matlab = MatlabCommand("ft_baseline_correction",
                               file_in,
                               tr,
                               cutoff_highpass)
        matlab.run()

        # read corrected time series
        arr = nb.load(os.path.join(path_file, "b{0}{1}".
                                   format(name_file, ext_file))).get_fdata()

    # temporal average
    arr_mean = np.mean(arr, axis=3) if mode == "mean" \
        else np.median(arr, axis=3)

    # write average
    data.header["dim"][0] = 3
    data.header["dim"][4] = 1
    output = nb.Nifti1Image(arr_mean, data.affine, data.header)
    file_out = os.path.join(path_output, "ssw.nii")
    nb.save(output, file_out)

    # bias field correction
    if apply_bias:
        n4 = N4BiasFieldCorrection()
        n4.inputs.dimension = 3
        n4.inputs.input_image = file_out
        n4.inputs.bias_image = os.path.join(path_output, "n4bias.nii")
        n4.inputs.output_image = os.path.join(path_output, "nssw.nii")
        n4.run()


if __name__ == "__main__":

    # description
    parser_description = "This program computes the temporal mean of an EPI " \
                         "time series. The time series is corrected for " \
                         "temporal drifts by applying a highpass filter and " \
                         "the temporal average is corrected for bias fields " \
                         "using ANTs. This average image can be used to " \
                         "show locations with pronounced static " \
                         "susceptibility effects (veins) in functional time " \
                         "series."

    input_help = "input file name of time series."
    tr_help = "volume repetition time (TR) in s."
    cutoff_help = "cutoff time in s for baseline correction (optional)."
    mode_help = "average mode (mean or median) (optional)."
    bias_help = "apply no bias field correction (optional)."

    # parse arguments from command line
    parser = argparse.ArgumentParser(description=parser_description)
    parser.add_argument('-i', '--in', type=str, help=input_help, dest='in_',
                        metavar="IN")
    parser.add_argument('-t', '--tr', type=float, help=tr_help)
    parser.add_argument('-c', '--cutoff', type=float, default=270,
                        help=cutoff_help)
    parser.add_argument('-m', '--mode', type=str, default="mean",
                        help=mode_help)
    parser.add_argument('-n', '--no_bias', action='store_false', help=bias_help)
    args = parser.parse_args()

    # run function
    static_susceptibility(args.in_, args.tr, args.cutoff, args.mode,
                          args.no_bias)
