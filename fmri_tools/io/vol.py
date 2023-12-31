# -*- coding: utf-8 -*-
"""Input/Output Nifti volume utilities."""

import subprocess

import nibabel as nb

from .filename import get_filename

__all__ = ["mri_convert", "copy_header"]


def mri_convert(file_in, file_out):
    """This function converts a volume file from freesurfer mgh/mgz to nifti nii/niigz
    format or vice versa.

    Parameters
    ----------
    file_in : str
        Full path of the input file.
    path_output : str
        Full path of the output file.

    Returns
    -------
    None.

    """
    # get filename
    _, _, ext_in = get_filename(file_in)
    _, _, ext_out = get_filename(file_out)

    if ext_in or ext_out not in [".nii", ".nii.gz", ".mgh", ".mgz"]:
        raise ValueError("Invalid file extension!")

    # convert volume to nifti format
    command = "mri_convert"
    command += " --in_type " + ext_in.replace(".", "")
    command += " --out_type " + ext_out.replace(".", "")
    command += " --input_volume " + file_in
    command += " --output_volume " + file_out

    print("Execute: " + command)
    try:
        subprocess.run([command], check=True)
    except subprocess.CalledProcessError:
        print("Execuation failed!")


def copy_header(file_in):
    """The function reads the header information from a nifti file and copies it into an
    empty header. This is done to omit unnecessary header extensions.

    Parameters
    ----------
    file_in : str
        filename of nifti volume.

    Returns
    -------
    header_new : niiheader
        Cleaned nibabel header.

    """
    # load header
    header = nb.load(file_in).header

    # initialize empty header
    header_new = nb.Nifti1Header()

    # copy header content
    header_new["sizeof_hdr"] = header["sizeof_hdr"]
    header_new["data_type"] = header["data_type"]
    header_new["db_name"] = header["db_name"]
    header_new["extents"] = header["extents"]
    header_new["session_error"] = header["session_error"]
    header_new["regular"] = header["regular"]
    header_new["dim_info"] = header["dim_info"]
    header_new["dim"] = header["dim"]
    header_new["intent_p1"] = header["intent_p1"]
    header_new["intent_p2"] = header["intent_p2"]
    header_new["intent_p3"] = header["intent_p3"]
    header_new["intent_code"] = header["intent_code"]
    header_new["datatype"] = header["datatype"]
    header_new["bitpix"] = header["bitpix"]
    header_new["slice_start"] = header["slice_start"]
    header_new["pixdim"] = header["pixdim"]
    header_new["vox_offset"] = header["vox_offset"]
    header_new["scl_slope"] = header["scl_slope"]
    header_new["scl_inter"] = header["scl_inter"]
    header_new["slice_end"] = header["slice_end"]
    header_new["slice_code"] = header["slice_code"]
    header_new["xyzt_units"] = header["xyzt_units"]
    header_new["cal_max"] = header["cal_max"]
    header_new["cal_min"] = header["cal_min"]
    header_new["slice_duration"] = header["slice_duration"]
    header_new["toffset"] = header["toffset"]
    header_new["glmax"] = header["glmax"]
    header_new["glmin"] = header["glmin"]
    header_new["descrip"] = header["descrip"]
    header_new["aux_file"] = header["aux_file"]
    header_new["qform_code"] = header["qform_code"]
    header_new["sform_code"] = header["sform_code"]
    header_new["quatern_b"] = header["quatern_b"]
    header_new["quatern_c"] = header["quatern_c"]
    header_new["quatern_d"] = header["quatern_d"]
    header_new["qoffset_x"] = header["qoffset_x"]
    header_new["qoffset_y"] = header["qoffset_y"]
    header_new["qoffset_z"] = header["qoffset_z"]
    header_new["srow_x"] = header["srow_x"]
    header_new["srow_y"] = header["srow_y"]
    header_new["srow_z"] = header["srow_z"]
    header_new["intent_name"] = header["intent_name"]
    header_new["magic"] = header["magic"]

    return header_new
