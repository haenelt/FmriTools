# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb


def get_weighted_vfs(input_vfs, input_snr, hemi, path_output):
    """Get weighted VFS.

    This function scales the visual fieldsign assignments using a given snr 
    estimate.    

    Parameters
    ----------
    input_vfs : str
        Input file with visual fieldsign assignments.
    input_snr : str
        Input file with snr estimates.
    hemi : str
        Hemisphere.
    path_output : str
        Output folder where output is saved.

    Returns
    -------
    None.
    
    """

    # load visual fieldsign map
    vfs_img = nb.load(input_vfs)
    vfs_array = vfs_img.get_fdata()

    # load snr map
    snr_img = nb.load(input_snr)
    snr_array = snr_img.get_fdata()

    # normalize snr map
    snr_array = snr_array / np.max(snr_array)

    # get weighted visual fieldsign map
    vfs_array = vfs_array * snr_array

    # write output
    out_img = nb.Nifti1Image(vfs_array, vfs_img.affine, vfs_img.header)
    nb.save(out_img, os.path.join(path_output, hemi + ".fieldsign_weighted.mgh"))
