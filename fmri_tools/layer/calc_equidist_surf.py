# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys

# external inputs
import numpy as np
import nibabel as nb
from nibabel.affines import apply_affine
from nibabel.freesurfer.io import read_geometry, write_geometry
from nighres.laminar import profile_meshing
from gbb.utils import vox2ras
from gbb.io import get_filename

# local inputs
from fmri_tools.surface import smooth_surface


def _get_meshlines(vp, vw):
    """
    Get mesh to visualize point-to-point connections between white and pial 
    surface.
    """
    
    # face of line
    f_new = [[0,1,0]]

    v_res = []
    f_res = []
    for i in range(len(vw)):
    
        # add vertices
        v_new = [list(vw[i]),list(vp[i])]
        
        # update resulting vertex and face list
        v_res.extend(v_new)
        f_res.extend(f_new)
    
        f_new += 2
    
    # vertices and faces as array
    v_res = np.array(v_res)
    f_res = np.array(f_res)
    
    return v_res, f_res
  

def calc_equidist_surf(input_mesh, input_boundaries, path_output, n_layer, 
                       n_smooth=2, n_crop=2):
    """ Calc equidist surf
    
    The function computes a set of matched surface meshs to create equidistant
    layers within the cerebral cortex. First, from a starting mesh and a 
    corresponding 4D levelset image, outer boundary surfaces are generated.  
    Optionally, the input levelset image can be cropped by n_crop first and last 
    image slices. This can be done if the levelset image shows some artifacts at 
    the edges. Furthermore, the computed outer boundary surfaces can be 
    smoothed. The final layers are then created by placing equidistant vertices 
    between both outer surfaces and saving them as separate mesh files.

    Parameters
    ----------
    input_mesh : str
        Filename of mesh.
    input_boundaries : str
        Filename of 4D levelset image.
    path_output : str
        Path where output is written.
    n_layer : int
        Number of layers.
    n_smooth : int, optional
        Number of smoothing iterations. The default is 2.
    n_crop : int, optional
        Number of omitted image slices. The default is 2.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 17-10-2020
    Last modified: 18-10-2020

    """
    
    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    
    # get hemisphere and basename of freesurfer surface
    _, hemi, name = get_filename(input_mesh)

    # check filename
    if not hemi[:2] == "lh" and not hemi[:2] == "rh":
        sys.exit("Could not identify hemi from filename!")
    
    # load data
    vtx, fac = read_geometry(input_mesh)
    level = nb.load(input_boundaries)
    
    # get ras <-> vox transformation
    vox2ras_tkr, ras2vox_tkr = vox2ras(input_boundaries)

    # transform vertices to voxel space
    vtx = apply_affine(ras2vox_tkr, vtx)
    
    # cut volume edges
    if n_crop:
        arr_level = level.get_fdata()
        arr_level = arr_level[n_crop:-n_crop,n_crop:-n_crop,n_crop:-n_crop,:]
        level = nb.Nifti1Image(arr_level, level.affine, level.header)
        vtx -= n_crop
    
    mesh = dict()
    mesh["points"] = vtx
    mesh["faces"] = fac
    
    # compute new wm and csf borders
    border = profile_meshing(level,
                             mesh,
                             save_data=False, 
                             overwrite=False, 
                             output_dir=None,
                             file_name=None)
    
    vtx_white = border["profile"][0]["points"]
    vtx_pial = border["profile"][-1]["points"]

    # resize vertices  
    if n_crop:
        vtx_pial += n_crop
        vtx_white += n_crop
    
    # transform vertices to ras space
    vtx_white = apply_affine(vox2ras_tkr, vtx_white)
    vtx_pial = apply_affine(vox2ras_tkr, vtx_pial)
    
    # write temporary output
    tmp = np.random.randint(0, 10, 5)
    tmp_string = ''.join(str(i) for i in tmp)
    
    tmp_white = os.path.join(path_output, "tmp_white_"+tmp_string)
    tmp_pial = os.path.join(path_output, "tmp_pial_"+tmp_string)
    
    write_geometry(tmp_white, vtx_white, fac)
    write_geometry(tmp_pial, vtx_pial, fac)
    
    # smooth output
    if n_smooth:
        smooth_surface(tmp_white, tmp_white, n_smooth)
        smooth_surface(tmp_pial, tmp_pial, n_smooth)
        
        vtx_white, fac_white = read_geometry(tmp_white)            
        vtx_pial, fac_pial = read_geometry(tmp_pial)
    
    # mesh lines
    vtx_lines, fac_lines = _get_meshlines(vtx_pial, vtx_white)
    write_geometry(os.path.join(path_output,"mesh_line"), vtx_lines, fac_lines)
    
    # get final layers
    for i in range(n_layer):
        vtx_layer = vtx_white + i/(n_layer-1) * (vtx_pial - vtx_white)
        filename_out = os.path.join(path_output, hemi+"."+name+"_layer_"+str(i))
        write_geometry(filename_out, vtx_layer, fac)    

    # clean
    os.remove(tmp_white)
    os.remove(tmp_pial)
