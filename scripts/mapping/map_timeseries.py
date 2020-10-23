# -*- coding: utf-8 -*-

# python standard library inputs
import os
import re
import glob
import multiprocessing

# external inputs
import numpy as np
import nibabel as nb
from nibabel.freesurfer.io import read_geometry
from joblib import Parallel, delayed

# local inputs
from fmri_tools.io import get_filename
from fmri_tools.io import write_hdf5
from fmri_tools.surface import deform_surface
from fmri_tools.mapping import map2surface


"""
Sample time series

The following script applies a transformation to a set of surface files and maps
timeseries data onto each surface. The resulting array vertex x timepoint x
surface (layer) is saved as hdf5 file. First, all surfaces which are found in 
the given list of folders `path_surf` are sorted. It is expected that all 
surfaces have the prefix lh or rh to indicate the hemisphere. Furthermore, it is 
expected that surfaces have a number in the basename and can therefore be 
numerically sorted by that number (e.g. numbers can indicate separate cortical 
layers). Each surface is transformed to the source space of the epi time series. 
Then, data from each time point is sampled. If a list of timeseries is given, 
each time series array will be saved in a separate hdf5 file.

created by Daniel Haenelt
Date created: 23-10-2020
Last modified: 23-10-2020
"""

# input
vol_in = ["/data/pt_01880/Experiment2_Rivalry/p1/odc/GE_EPI1/Run_1/uadata.nii",
          "/data/pt_01880/Experiment2_Rivalry/p1/odc/GE_EPI1/Run_2/uadata.nii",
          ]
path_surf = ["/data/pt_01880/Experiment2_Rivalry/p1/anatomy/surf/lh/layer",
             "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/surf/rh/layer",
             ]
source2target_in = "/data/pt_01880/Experiment2_Rivalry/p1/deformation/odc/GE_EPI1/source2target.nii.gz"
interp_method = "trilinear"

# do not edit below

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    """
    see:
    https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-
    with-a-number-inside
    """
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def do_mapping(i, file_vol, file_surf, path_output, interp_method):
    """
    map on surface    
    """

    # get volume information
    aff = nb.load(file_vol).affine
    head = nb.load(file_vol).header
    head["dim"][0] = 3
    head["dim"][4] = 1

    # temporary filename
    tmp = np.random.randint(0, 10, 5)
    tmp_string = ''.join(str(x) for x in tmp)
    file_tmp = os.path.join(path_output, "tmp_"+tmp_string+".nii")

    output = nb.Nifti1Image(nb.load(file_vol).get_fdata()[:,:,:,i], aff, head)
    nb.save(output, file_tmp)   

    # do mapping                
    arr_map, aff_map, head_map = map2surface(input_surf=file_surf, 
                                             input_vol=file_tmp, 
                                             write_output=False,
                                             path_output="", 
                                             interp_method=interp_method,
                                             input_surf_target=None, 
                                             input_ind=None, 
                                             cleanup=True)  
    
    # remove temporary volume
    os.remove(file_tmp)
    
    return arr_map, aff_map, head_map

# number of cores
num_cores = multiprocessing.cpu_count()

# sort surface filenames
tmp = []
for i in range(len(path_surf)):
    tmp.extend(glob.glob(os.path.join(path_surf[i],"*")))

surf_left = []
surf_right = []
for i in tmp[:]:
    basename = os.path.basename(i)
    if basename.startswith("lh"):
        surf_left.append(i)
    elif basename.startswith("rh"):
        surf_right.append(i)

surf_left.sort(key=natural_keys)
surf_right.sort(key=natural_keys)

file_surf = []
file_surf.append(surf_left)
file_surf.append(surf_right)

# start mapping
for i in range(len(vol_in)):
    
    # get filename
    path_vol, name_vol, _ = get_filename(vol_in[i])
   
    # make output folder
    path_output = os.path.join(path_vol, "sampled")
    if not os.path.exists(path_output):
        os.makedirs(path_output)    

    # get volume information
    n_time = nb.load(vol_in[i]).header["dim"][4]
    
    for j in range(len(file_surf)):

        vtx, _ = read_geometry(file_surf[j][0])
        n_vtx = len(vtx)
        n_layer = len(file_surf[j])
        arr_res = np.zeros((n_vtx, n_time, n_layer))
        for l in range(n_layer):
            
            # deform mesh
            deform_surface(input_surf=file_surf[j][l],
                           input_orig=source2target_in,
                           input_deform=source2target_in,
                           input_target=vol_in[i],
                           path_output=path_output,
                           input_mask=None,
                           interp_method="trilinear",
                           smooth_iter=0,
                           flip_faces=False,
                           cleanup=True)

            # temporary surface
            file_def = os.path.join(path_output, 
                                    os.path.basename(file_surf[j][l])+"_def")    

            tmp = Parallel(n_jobs=num_cores)(
                delayed(do_mapping)(
                    t,
                    vol_in[i],
                    file_def, 
                    path_output, 
                    interp_method) for t in range(n_time)
                )

            for t in range(n_time):
                arr_res[:,t,l] = tmp[t][0]
                affine = tmp[t][1]
                header = tmp[t][2]
            
            # remove deformed surface
            os.remove(file_def)
        
        # write hdf5
        _, hemi, _ = get_filename(file_surf[j][0])
        file_out = os.path.join(path_output, 
                                hemi+"."+name_vol+"_nlayer"+str(n_layer)+".hdf5")
        write_hdf5(file_out, arr_res, affine, header)
