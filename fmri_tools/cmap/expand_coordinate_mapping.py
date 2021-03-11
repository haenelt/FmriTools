# -*- coding: utf-8 -*-

# python standard library inputs
import os
import random

# external inputs
import numpy as np
import nibabel as nb

# local inputs
from fmri_tools.io.get_filename import get_filename
from fmri_tools.cmap.generate_coordinate_mapping import generate_coordinate_mapping
from fmri_tools.utils.apply_affine_chunked import apply_affine_chunked


def expand_coordinate_mapping(cmap_in, path_output=None, name_output=None, 
                              write_output=False):
    """ Expand coordinate mapping
    
    This function removes black background in a coordinate mapping to omit 
    interpolation problems at the edges of a coordinate slab within a larger 
    volume. Based on the cmap, a transformation matrix is computed from randomly 
    sampled data points within the slab. The transformation matrix is then 
    applied to all background voxels. Hence, this method is only really precise 
    for coordinate mappings representing an affine transformation. However, this 
    function can also be applied to nonlinear coordinate mappings since the 
    preliminary goal is to avoid problems at the slab edges. Therefore, the 
    actual data sampling should not be affected. The code snippet for computing 
    the transformation matrix is taken from [1].    
    

    Parameters
    ----------
    cmap_in : str
        Filename of coordinate mapping.
    path_output : str, optional
        Path where output is written. The default is None.
    name_output : str, optional
        Basename of output volume. The default is None.
    write_output : bool, optional
        Write nifti volume. The default is False.

    Returns
    -------
    output : niimg
        Corrected coordinate mapping.

    References
    -------
    .. [1] https://stackoverflow.com/questions/56220626/how-to-compute-
    conformal-affine-transformation
    
    Notes
    -------
    created by Daniel Haenelt
    Date created: 18-06-2020
    Last modified: 11-03-2021
    
    """
    
    # get file extension of cmap
    _, _, ext_cmap = get_filename(cmap_in)
    
    # load target cmap and generate source cmap 
    cmap_target = nb.load(cmap_in)
    cmap_source = generate_coordinate_mapping(cmap_in, pad=0)
    
    arr_cmap_target = cmap_target.get_fdata()
    arr_cmap_source = cmap_source.get_fdata()
    
    # get image dimensions
    xdim = cmap_source.header["dim"][1]
    ydim = cmap_source.header["dim"][2]
    zdim = cmap_source.header["dim"][3]
    
    # random selection of 4 data points
    s_coords = []
    t_coords = []
    
    pts = np.where(arr_cmap_target[:,:,:,0] != 0)
    r = random.sample(np.arange(len(pts[0])).tolist(),len(pts[0]))[:4]
    
    s_coords.append([pts[0][r[0]], pts[1][r[0]], pts[2][r[0]]])
    s_coords.append([pts[0][r[1]], pts[1][r[1]], pts[2][r[1]]])
    s_coords.append([pts[0][r[2]], pts[1][r[2]], pts[2][r[2]]])
    s_coords.append([pts[0][r[3]], pts[1][r[3]], pts[2][r[3]]])
    
    t_coords.append(arr_cmap_target[s_coords[0][0], s_coords[0][1], s_coords[0][2],:].tolist())
    t_coords.append(arr_cmap_target[s_coords[1][0], s_coords[1][1], s_coords[1][2],:].tolist())
    t_coords.append(arr_cmap_target[s_coords[2][0], s_coords[2][1], s_coords[2][2],:].tolist())
    t_coords.append(arr_cmap_target[s_coords[3][0], s_coords[3][1], s_coords[3][2],:].tolist())
    
    # get transformation matrix (target -> source)
    l = len(t_coords)
    B = np.vstack([np.transpose(t_coords), np.ones(l)])
    D = 1.0 / np.linalg.det(B)
    
    entry = lambda r,d: np.linalg.det(np.delete(np.vstack([r, B]), (d+1), axis=0))
    M = [[(-1)**i * D * entry(R, i) for i in range(l)] for R in np.transpose(s_coords)]
    A, t = np.hsplit(np.array(M), [l-1])
    t = np.transpose(t)[0]
    
    # unittests
    print("Test cmap expansion:")
    for p, P in zip(np.array(t_coords), np.array(s_coords)):
      image_p = np.dot(A, p) + t
      result = "[OK]" if np.allclose(image_p, P) else "[ERROR]"
      print(p, " mapped to: ", image_p, " ; expected: ", P, result)
      
    # get affine transformation matrix by adding translation vector
    M = np.zeros((4,4))
    M[:3,:3] = A
    M[:3,-1] = t
    M[-1,-1] = 1
    
    # get final transformation matrix (source -> target)
    M = np.linalg.inv(M)
    
    # transform source volume
    x = arr_cmap_source[:,:,:,0].flatten()
    y = arr_cmap_source[:,:,:,1].flatten()
    z = arr_cmap_source[:,:,:,2].flatten()
    
    source_listed = np.array([x,y,z]).T
    source_transformed = apply_affine_chunked(M, source_listed)
    
    x_new = np.reshape(source_transformed[:,0], (xdim,ydim,zdim))
    y_new = np.reshape(source_transformed[:,1], (xdim,ydim,zdim))
    z_new = np.reshape(source_transformed[:,2], (xdim,ydim,zdim))

    # overwrite new cmap with old coordinate mapping (so that only background remains)    
    x_new[arr_cmap_target[:,:,:,0] > 0] = arr_cmap_target[:,:,:,0][arr_cmap_target[:,:,:,0] > 0]
    y_new[arr_cmap_target[:,:,:,1] > 0] = arr_cmap_target[:,:,:,1][arr_cmap_target[:,:,:,1] > 0]
    z_new[arr_cmap_target[:,:,:,2] > 0] = arr_cmap_target[:,:,:,2][arr_cmap_target[:,:,:,2] > 0]
    
    # overwrite input cmap array with final cmap array
    arr_cmap_target[:,:,:,0] = x_new
    arr_cmap_target[:,:,:,1] = y_new
    arr_cmap_target[:,:,:,2] = z_new
    
    # nibabel instance of final cmap
    output = nb.Nifti1Image(arr_cmap_target, cmap_target.affine, cmap_target.header)
    
    # write output
    if write_output:
        nb.save(output, os.path.join(path_output, name_output+ext_cmap))
    
    return output
