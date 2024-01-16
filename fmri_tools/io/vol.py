# -*- coding: utf-8 -*-
"""Input/Output Nifti volume utilities."""

import os
import subprocess

import nibabel as nb
import numpy as np

from .filename import get_filename

__all__ = ["load_volume", "save_volume", "mri_convert", "copy_header", "surface_voxel"]


def load_volume(volume):
    """Load volumetric data into a
    `Nibabel SpatialImage <http://nipy.org/nibabel/reference/nibabel.spatialimages.html#nibabel.spatialimages.SpatialImage>`_

    Parameters
    ----------
    volume: niimg
        Volumetric data to be loaded, can be a path to a file that nibabel can
        load, or a Nibabel SpatialImage

    Returns
    ----------
    image: Nibabel SpatialImage

    Notes
    ----------
    Originally created as part of Laminar Python [1]_ .

    References
    -----------
    .. [1] Huntenburg et al. (2017), Laminar Python: Tools for cortical
       depth-resolved analysis of high-resolution brain imaging data in
       Python. DOI: 10.3897/rio.3.e12346
    """
    # if input is a filename, try to load it
    # python 2 version if isinstance(volume, basestring):
    if isinstance(volume, str):
        # importing nifti files
        image = nb.load(volume)
    # if volume is already a nibabel object
    elif isinstance(volume, nb.spatialimages.SpatialImage):
        image = volume
    else:
        raise ValueError(
            "Input volume must be a either a path to a file in a "
            "format that Nibabel can load, or a nibabel"
            "SpatialImage."
        )
    return image


def save_volume(filename, volume, dtype="float32", overwrite_file=True):
    """Save volumetric data that is a
    `Nibabel SpatialImage <http://nipy.org/nibabel/reference/nibabel.spatialimages.html#nibabel.spatialimages.SpatialImage>`_
    to a file

    Parameters
    ----------
    filename: str
        Full path and filename under which volume should be saved. The
        extension determines the file format (must be supported by Nibabel)
    volume: Nibabel SpatialImage
        Volumetric data to be saved
    dtype: str, optional
        Datatype in which volumetric data should be stored (default is float32)
    overwrite_file: bool, optional
        Overwrite existing files (default is True)

    Notes
    ----------
    Originally created as part of Laminar Python [1]_ .

    References
    -----------
    .. [1] Huntenburg et al. (2017), Laminar Python: Tools for cortical
       depth-resolved analysis of high-resolution brain imaging data in
       Python. DOI: 10.3897/rio.3.e12346
    """  # noqa
    if dtype is not None:
        volume.set_data_dtype(dtype)
    if os.path.isfile(filename) and overwrite_file is False:
        print(
            "\nThis file exists and overwrite_file was set to False, " "file not saved."
        )
    else:
        try:
            volume.to_filename(filename)
            print("\nSaving {0}".format(filename))
        except AttributeError:
            print("\nInput volume must be a Nibabel SpatialImage.")


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


def surface_voxel(file_in, path_output):
    """Compute surface voxels for to illustrate curvature dependence in surface
    representations as introduced in [1].

    Parameters
    ----------
    file_in : str
        Input image.
    path_output : str
        Path where to write output image.

    Returns
    -------
    None.

    References
    -------
    .. [1] Kay, K, et al. A critical assessment of data quality and venous
    effects in ultra-high-resolution fMRI, bioRxiv, 1--45 (2018).

    """
    # parameters
    x_calc = 1
    y_calc = 1
    z_calc = 1
    x_max = 1
    y_max = 2
    z_max = 4

    # make subfolders
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # load img
    img = nb.load(file_in)

    # img range
    x_size = img.shape[0]
    y_size = img.shape[1]
    z_size = img.shape[2]

    # compute surface voxels
    img_res = np.zeros([x_size, y_size, z_size])
    for i in range(x_size):
        for j in range(y_size):
            for k in range(z_size):
                if np.mod(i, 2) != 0 and x_calc != 0:
                    img_res[i, j, k] += x_max
                if np.mod(j, 2) != 0 and y_calc != 0:
                    img_res[i, j, k] += y_max
                if np.mod(k, 2) != 0 and z_calc != 0:
                    img_res[i, j, k] += z_max

    # write output image
    newimg = nb.Nifti1Image(img_res, img.affine, img.header)
    nb.save(newimg, os.path.join(path_output, "surface_voxel.nii"))
