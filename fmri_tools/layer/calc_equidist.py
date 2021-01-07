# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys
from collections import Counter

# external inputs
import numpy as np
import nibabel as nb
from nibabel.affines import apply_affine
from nibabel.freesurfer.io import read_geometry
from skimage import measure
from nighres.surface import probability_to_levelset
from scipy.ndimage.morphology import binary_fill_holes
from gbb.utils.vox2ras import vox2ras

# local inputs
from fmri_tools.io import get_filename
from fmri_tools.utils import resample_volume
from fmri_tools.surface import upsample_surf_mesh


def calc_equidist(input_white, input_pial, input_vol, n_layers, path_output,
                  r=[0.4,0.4,0.4], n_iter=2, pathLAYNII="", debug=False):    
    """ Calc equidist
    
    This function computes equidistant layers in volume space from input pial 
    and white surfaces in freesurfer format using the laynii function 
    LN_GROW_LAYERS. The input surfaces do not have to cover the whole brain. 
    Number of vertices and indices do not have to correspond between surfaces.    

    Parameters
    ----------
    input_white : str
        Filename of white surface.
    input_pial : str
        Filename of pial surface.
    input_vol : str
        Filename of reference volume.
    n_layers : int
        Number of generated layers.
    path_output : str
        Path where output is written.
    r : list, optional
        Array of new voxel sizes for reference volume upsampling (if not None). 
        The default is [0.4,0.4,0.4].
    n_iter : int, optional
        Number of surface upsampling iterations. The default is 2.    
    pathLAYNII : str, optional
        Path to laynii folder. The default is "".
    debug : bool, optional
        Write out some intermediate files. The default is False. 

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 31-05-2020
    Last modified: 17-10-2020   

    """
    
    # add laynii to search path
    if len(pathLAYNII):
        os.environ["PATH"] += os.path.pathsep + pathLAYNII
    
    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    
    # get hemi from filename
    _, hemi, _ = get_filename(input_white)
    if not hemi[:2] == "lh" and not hemi[:2] == "rh":
        sys.exit("Could not identify hemi from filename!")
    
    # new filenames in output folder
    res_white = os.path.join(path_output,hemi+".white")
    res_pial = os.path.join(path_output,hemi+".pial")
    res_vol = os.path.join(path_output,"epi_upsampled.nii")
    
    # upsample reference volume and input surface
    resample_volume(input_vol, res_vol, dxyz=r, rmode="Cu")    
    upsample_surf_mesh(input_white, res_white, n_iter, "linear")
    upsample_surf_mesh(input_pial, res_pial, n_iter, "linear")
    
    # get affine ras2vox-tkr transformation to reference volume
    _, ras2vox_tkr = vox2ras(res_vol)
    
    # load surface
    vtx_white, _ = read_geometry(res_white) 
    vtx_pial, _ = read_geometry(res_pial)
    
    # load volume
    vol = nb.load(res_vol)
    
    # apply ras2vox to coords
    vtx_white = np.round(apply_affine(ras2vox_tkr, vtx_white)).astype(int)
    vtx_pial = np.round(apply_affine(ras2vox_tkr, vtx_pial)).astype(int)
    
    # surfaces to lines in volume
    white_array = np.zeros(vol.header["dim"][1:4])
    white_array[vtx_white[:,0],vtx_white[:,1],vtx_white[:,2]] = 1   
    white = nb.Nifti1Image(white_array, vol.affine, vol.header)
    
    pial_array = np.zeros(vol.header["dim"][1:4])
    pial_array[vtx_pial[:,0],vtx_pial[:,1],vtx_pial[:,2]] = 1
    pial = nb.Nifti1Image(pial_array, vol.affine, vol.header)
    
    # make wm
    white_label_array = np.zeros_like(white_array)
    for i in range(np.shape(white_label_array)[2]):
        white_label_array[:,:,i] = binary_fill_holes(white_array[:,:,i])
    white_label_array = white_label_array - white_array
    white_label_array = measure.label(white_label_array, connectivity=1)
    white_label_flatten = np.ndarray.flatten(white_label_array)
    white_label_flatten = white_label_flatten[white_label_flatten > 0]
    label_number = Counter(white_label_flatten).most_common(1)[0][0]
    white_label_array[white_label_array != label_number] = 0
    white_label_array[white_label_array > 0] = 1    
    white_label = nb.Nifti1Image(white_label_array, vol.affine, vol.header)
    
    # make csf
    pial_label_array = np.zeros_like(pial_array)
    for i in range(np.shape(pial_label_array)[2]):
        pial_label_array[:,:,i] = binary_fill_holes(pial_array[:,:,i])
    pial_label_array = pial_label_array - pial_array
    pial_label_array = measure.label(pial_label_array, connectivity=1)
    pial_label_flatten = np.ndarray.flatten(pial_label_array)
    pial_label_flatten = pial_label_flatten[pial_label_flatten > 0]
    label_number = Counter(pial_label_flatten).most_common(1)[0][0]
    pial_label_array[pial_label_array != label_number] = 0
    pial_label_array[pial_label_array > 0] = 1    
    pial_label_array = pial_label_array + pial_array # add csf line again
    pial_label = nb.Nifti1Image(pial_label_array, vol.affine, vol.header)
       
    # make gm
    ribbon_label_array = pial_label_array - white_label_array
    ribbon_label_array[ribbon_label_array != 1] = 0
    ribbon_label = nb.Nifti1Image(ribbon_label_array, vol.affine, vol.header)
    
    # make rim
    rim_array = np.zeros_like(ribbon_label_array)
    rim_array[ribbon_label_array == 1] = 3
    rim_array[pial_array == 1] = 1
    rim_array[white_array == 1] = 2
       
    output = nb.Nifti1Image(rim_array, vol.affine, vol.header)
    nb.save(output, os.path.join(path_output, "rim.nii"))    
    
    # grow layers using laynii
    vinc = 40
    os.system("LN_GROW_LAYERS" + \
              " -rim " + os.path.join(path_output, "rim.nii") + \
              " -vinc " + str(vinc) + \
              " -N " + str(n_layers) + \
              " -threeD" + \
              " -output " + os.path.join(path_output, "layers.nii"))

    # tranform label to levelset    
    binary_array = white_label_array + ribbon_label_array + pial_label_array
    binary_array[binary_array != 0] = 1 
    
    layer_array =  nb.load(os.path.join(path_output, "layers.nii")).get_fdata()
    layer_array += 1
    layer_array[layer_array == 1] = 0
    layer_array[white_label_array == 1] = 1 # fill wm
    
    if debug:
        out_debug = nb.Nifti1Image(layer_array, vol.affine, vol.header)
        nb.save(out_debug, 
                os.path.join(path_output, "layers_plus_white_debug.nii"))
    
    level_array = np.zeros(np.append(vol.header["dim"][1:4],n_layers))
    for i in range(n_layers):
        print("Probabilty to levelset for layer: "+str(i+1))
        
        temp_layer_array = binary_array.copy()
        temp_layer_array[layer_array > i+1] = 0
        temp_layer = nb.Nifti1Image(temp_layer_array, vol.affine, vol.header)
    
        # write control output
        if debug:
            nb.save(temp_layer, 
                    os.path.join(path_output, "layer_"+str(i)+"_debug.nii"))
    
        # transform binary image to levelset image
        res = probability_to_levelset(temp_layer)
        
        # sort levelset image into 4d array
        level_array[:,:,:,i] = res["result"].get_fdata()

    # move border to voxel center
    level_array -= 0.5

    # levelset image
    vol.header["dim"][0] = 4
    vol.header["dim"][4] = n_layers
    levelset = nb.Nifti1Image(level_array, vol.affine, vol.header)
    
    # write niftis
    nb.save(white, os.path.join(path_output, "wm_line.nii"))
    nb.save(pial, os.path.join(path_output, "csf_line.nii"))
    nb.save(white_label, os.path.join(path_output, "wm_label.nii"))
    nb.save(pial_label, os.path.join(path_output, "csf_label.nii"))
    nb.save(ribbon_label, os.path.join(path_output, "gm_label.nii"))
    nb.save(levelset, os.path.join(path_output, "boundaries.nii"))    
