# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb
from nibabel.freesurfer.io import read_geometry, write_morph_data
from gbb.normal import get_normal
from gbb.utils.vox2ras import vox2ras

# local inputs
from ..io.get_filename import get_filename
from ..utils.apply_affine_chunked import apply_affine_chunked


def b0_orientation(surf_in, vol_in, write_output=False, path_output="",
                   name_output=""):
    """B0 orientation.
    
    This function computes the angle between surface normals and B0-direction 
    per vertex.    

    Parameters
    ----------
    surf_in : str
        Input of surface mesh.
    vol_in : str
        Input of corresponding nifti volume.
    write_output : bool, optional
        Write out to disk (boolean). The default is False.
    path_output : str, optional
        Path where to save output. The default is "".
    name_output : str, optional
        Basename of output file. The default is "".

    Returns
    -------
    theta : ndarray
        Angle in radians.
    
    """
    
    # make subfolders
    if write_output and not os.path.exists(path_output):
        os.makedirs(path_output)
    
    # get hemi from surface filename
    _, hemi, _ = get_filename(surf_in)
    
    # load surface
    vtx, fac = read_geometry(surf_in)
    
    # get transformation matrix
    _, r2v = vox2ras(vol_in)      # ras-tkr -> voxel
    v2s = nb.load(vol_in).affine  # voxel -> scanner-ras
    m = v2s.dot(r2v)
    
    # apply affine transformation
    vtx = apply_affine_chunked(m, vtx)
    
    # get surface normals
    n = get_normal(vtx, fac)
    
    # get angle between b0 and surface normals in radians    
    theta = np.arccos(np.dot(n, [0, 0, 1]))
    
    # write output
    if write_output:
        write_morph_data(os.path.join(path_output, hemi+"."+name_output), theta)

    return theta
