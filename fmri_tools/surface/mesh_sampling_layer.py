# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys
import shutil as sh
from os.path import join, exists, basename, splitext

# external inputs
import numpy as np
import nibabel as nb
from nighres.laminar import profile_sampling

# local inputs
from fmri_tools.io.get_filename import get_filename
from fmri_tools.utils.upsample_volume import upsample_volume
from fmri_tools.mapping import map2surface


def mesh_sampling_layer(surf_in, file_in, boundaries_in, path_output, layer, 
                        r=[0.4,0.4,0.4], interpolation="Cu", 
                        average_layer=False, write_profile=True, 
                        write_upsampled=True):
    """ Mesh sampling layer

    This function samples data from an image volume to a surface mesh from 
    specific layers defined by a levelset image. If average_layer is true, the 
    parameter layer should contain only two integers which denote the start and 
    ending layer.    

    Parameters
    ----------
    surf_in : str
        Filename of input surface mesh.
    file_in : str
        Filename of input volume from which data is sampled.
    boundaries_in : str
        Filename of 4D levelset image.
    path_output : str
        Path where output is written.
    layer : list
        Which layers to sample (array of integers).
    r : list, optional
        Destination voxel size after upsampling (performed if not None). The 
        default is [0.4,0.4,0.4].
    interpolation : str, optional
        Interpolation method for upsampling of file from which data is sampled. 
        The default is "Cu".
    average_layer : bool, optional
        Average across cortex. The default is False.
    write_profile : bool, optional
        Write sampled profile. The default is True.
    write_upsampled : bool, optional
        Write upsampled file. The default is True.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 18-12-2019
    Last modified: 19-10-2020

    """
    
    # make output folder
    if not exists(path_output):
        os.makedirs(path_output)
    
    # filenames
    _, name_file, ext_file = get_filename(file_in)   
    _, hemi, name_surf = get_filename(surf_in)  
     
    name_surf = name_surf[1:]
    name_profile = splitext(basename(file_in))[0]+"_profile"
    
    # check hemi
    if not hemi == "lh" and not hemi == "rh":
        sys.exit("Could not identify hemi from filename!")
    
    # upsample volume
    if not r == None:
        name_file = name_file+"_upsampled"
        upsample_volume(file_in, join(path_output, name_file+ext_file), r, interpolation)
    else:
        if file_in != join(path_output, name_file+ext_file):
            sh.copyfile(file_in, join(path_output, name_file+ext_file))
        
    # get profile sampling            
    tmp = np.random.randint(0, 10, 5)
    tmp_string = ''.join(str(i) for i in tmp)
    profile = profile_sampling(boundaries_in, 
                               join(path_output, name_file+ext_file),
                               save_data=write_profile, 
                               overwrite=write_profile,
                               output_dir=path_output,
                               file_name="profile_"+tmp_string)
        
    # rename profile sampling output
    if write_profile:
        os.rename(join(path_output, "profile_"+tmp_string+"_lps-data.nii.gz"),
                  join(path_output, name_profile+".nii.gz"))
    
    # load profile
    if write_profile:
        data = nb.load(join(path_output, name_profile+".nii.gz"))
    else:
        data = profile["result"]
    data.header["dim"][0] = 3        
    
    # do mapping
    tmp2 = np.random.randint(0, 10, 5)
    tmp2_string = ''.join(str(i) for i in tmp2)
    if not average_layer:
        
        for i in range(len(layer)):
            data_array = data.get_fdata()[:,:,:,layer[i]]
            out = nb.Nifti1Image(data_array, data.affine, data.header)
            nb.save(out, join(path_output,"temp_"+tmp2_string+".nii"))
            
            # do the mapping
            map2surface(surf_in, 
                        join(path_output,"temp_"+tmp2_string+".nii"),
                        path_output,
                        interp_method="nearest",
                        input_white=None, 
                        input_ind=None, 
                        cleanup=True)

            # rename mapping file
            os.rename(join(path_output,hemi+".temp_"+tmp2_string+"_"+name_surf+".mgh"),
                      join(path_output,hemi+"."+name_file+"_layer"+str(layer[i])+".mgh"))

    else:

        if len(layer) != 2:
            sys.exit("For averaging, layer should only contain two elements!")
        
        data_array = data.get_fdata()[:,:,:,layer[0]:layer[1]]
        data_array = np.mean(data_array, axis=3)
        out = nb.Nifti1Image(data_array, data.affine, data.header)
        nb.save(out, join(path_output,"temp_"+tmp2_string+".nii"))
                
        # do the mapping
        map2surface(surf_in, 
                    join(path_output,"temp_"+tmp2_string+".nii"),
                    path_output,
                    interp_method="nearest",
                    input_white=None, 
                    input_ind=None, 
                    cleanup=True)

        # rename mapping file
        os.rename(join(path_output,hemi+".temp_"+tmp2_string+"_"+name_surf+".mgh"),
                  join(path_output,hemi+"."+name_file+"_avg_layer"+str(layer[0])+"_"+str(layer[1])+".mgh"))
        
    # clean temp
    os.remove(join(path_output,"temp_"+tmp2_string+".nii"))
    
    # clean file
    if not write_upsampled:
        os.remove(join(path_output, name_file+ext_file))
