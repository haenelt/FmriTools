# -*- coding: utf-8 -*-

# python standard library inputs
import os
import shutil as sh

# external inputs
import numpy as np
import nibabel as nb
from nibabel.affines import apply_affine
from nighres.laminar import profile_sampling
from gbb.utils.vox2ras import vox2ras

# local inputs
from fmri_tools.cmap.generate_coordinate_mapping import generate_coordinate_mapping
from fmri_tools.utils.upsample_volume import upsample_volume


def get_thickness(boundaries_in, ref_in, hemi, path_output, r=[0.4,0.4,0.4]):
    """ Get thickness

    This function computes the cortical thickness as euclidean distances between 
    vertex ras coordinates from outer levelset boundaries.    

    Parameters
    ----------
    boundaries_in : str
        Filename of 4D levelset image.
    ref_in : str
        Filename of reference volume for coordinate transformation.
    hemi : str
        Hemisphere.
    path_output : str
        Path where output is written.
    r : list, optional
        Destination voxel size after upsampling (performed if not None). The 
        default is [0.4,0.4,0.4].

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 18-12-2019
    Last modified: 13-10-2020
    
    """
    
    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    
    # upsample volume
    if not r == None:
        upsample_volume(ref_in, os.path.join(path_output, "ref.nii"), r, "Cu")
    else:
        sh.copyfile(ref_in, os.path.join(path_output, "ref.nii"))
    
    # make coordinate mapping
    cmap = generate_coordinate_mapping(boundaries_in, pad=0)
    cmap.header["dim"][0] = 1
    
    # get voxel to vertex ras coordinate transformation
    vox2ras_tkr, _ = vox2ras(os.path.join(path_output, "ref.nii"))
    
    # apply transformation to cmap
    ras_array = apply_affine(vox2ras_tkr, cmap.get_fdata())
    
    # split coordinates into single dimensions
    x_ras = nb.Nifti1Image(ras_array[:,:,:,0], cmap.affine, cmap.header)
    y_ras = nb.Nifti1Image(ras_array[:,:,:,1], cmap.affine, cmap.header)
    z_ras = nb.Nifti1Image(ras_array[:,:,:,2], cmap.affine, cmap.header)
    
    # get profile sampling
    x_profile = profile_sampling(boundaries_in, x_ras)
    y_profile = profile_sampling(boundaries_in, y_ras)
    z_profile = profile_sampling(boundaries_in, z_ras)
        
    # compute euclidean distance between outer levelset boundaries
    x_array = x_profile["result"].get_fdata()
    y_array = y_profile["result"].get_fdata()
    z_array = z_profile["result"].get_fdata()
    
    x_diff_array = np.square( x_array[:,:,:,-1] - x_array[:,:,:,0] )
    y_diff_array = np.square( y_array[:,:,:,-1] - y_array[:,:,:,0] )
    z_diff_array = np.square( z_array[:,:,:,-1] - z_array[:,:,:,0] )
    
    r = np.sqrt( x_diff_array + y_diff_array + z_diff_array )
    
    # set unrealistic values to zero
    r[r>10] = 0
    
    # hemi suffix
    if hemi == "lh":
        hemi_suffix = "_left"
    elif hemi == "rh":
        hemi_suffix = "_right"
    else:
        hemi_suffix = None
    
    # write nifti
    output = nb.Nifti1Image(r, cmap.affine, cmap.header)
    nb.save(output, os.path.join(path_output, "thickness"+hemi_suffix+".nii"))
    
    # remove temporary reference volume
    os.remove(os.path.join(path_output, "ref.nii"))
