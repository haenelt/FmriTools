# -*- coding: utf-8 -*-

# python standard library inputs
import argparse

# external inputs
import numpy as np
import nibabel as nb


def calc_mip(file_in, file_out, size, axis=2, mode="min"):
    """Minimum/Maximum intensity projection of a 3D nifti image. For each slice along
    the specified axis, the maximum or minimum intensity is computed from a number of
    neighboring slices around the current slice. The number of neighboring slices in
    each direction is defined by the size parameter.

    Parameters
    ----------
    file_in : str
        File name of input image.
    file_out : str
        File name of output image.
    size : int
        Number of included slices in each direction.
    axis : int, optional
        Axis along which to compute the intensity projection.
    mode : str, optional
        Mode (min, max).

    Returns
    -------
    None.

    """

    if axis < 0 or axis > 2:
        raise ValueError("Axis must be between 0 and 2.")

    # load input image
    data = nb.load(file_in)
    dim = data.header["dim"][axis + 1]
    arr = data.get_fdata()

    arr = np.moveaxis(arr, axis, 0)
    res = np.zeros_like(arr)
    for i in range(dim):
        index_low = i - size
        if index_low < 0:
            index_low = 0

        index_high = i + size
        if index_high > dim:
            index_high = dim

        if mode == "min":
            res[i,:,:] = np.min(arr[index_low:index_high,:,:], axis=0)
        elif mode == "max":
            res[i,:,:] = np.max(arr[index_low:index_high,:,:], axis=0)
        else:
            raise ValueError("Mode not supported")

    res = np.moveaxis(res, 0, axis)
    output = nb.Nifti1Image(res, data.affine, data.header)
    nb.save(output, file_out)


if __name__ == "__main__":

    # description
    parser_description = "This program computes a minimum/maximum intensity " \
                         "projection of a 3D nifti image."

    in_help = "input file name."
    out_help = "output file name."
    size_help = "number of included slices in each direction."
    axis_help = "projection axis."
    mode_help = "mode (min, max)."

    # parse arguments from command line
    parser = argparse.ArgumentParser(description=parser_description)
    parser.add_argument('-i', '--in', type=str, help=in_help, dest='in_', metavar="IN")
    parser.add_argument('-o', '--out', type=str, help=out_help)
    parser.add_argument('-s', '--size', type=int, help=size_help)
    parser.add_argument('-a', '--axis', type=int, help=axis_help)
    parser.add_argument('-m', '--mode', type=str, help=mode_help)
    args = parser.parse_args()

    # run function
    calc_mip(args.in_, args.out, args.size, args.axis, args.mode)
