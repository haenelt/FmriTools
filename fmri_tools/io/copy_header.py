# -*- coding: utf-8 -*-

# external inputs
import nibabel as nb


def copy_header(file_in):
    """
    The function reads the header information from a nifti file and copies it 
    into an empty header. This is done to omit unnecessary header extensions.
    Inputs:
        *file_in: filename of nifti volume.
    Outputs:
        *header_new: cleaned nibabel header.

    created by Daniel Haenelt
    Date created: 23-06-2020             
    Last modified: 12-10-2020
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
