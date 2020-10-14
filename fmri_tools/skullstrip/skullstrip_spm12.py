# -*- coding: utf-8 -*-

# python standard library inputs
import os
import shutil

# external inputs
import nibabel as nb
from nipype.interfaces.spm import NewSegment
from nipype.interfaces.matlab import MatlabCommand


def skullstrip_spm12(filename, pathSPM12, path_output):
    """ Skullstrip SPM12

    The computation of the skullstrip mask is done on the PD-weighted INV2 
    image. According to S. Kashyap, this shows the best results. Outputs are 
    written in a subfolder of the given output path.    

    Parameters
    ----------
    filename : str
        Path of input image.
    pathSPM12 : str
        Path to SPM12 toolbox.
    path_output : str
        Path where output is saved.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 01-11-2018             
    Last modified: 12-10-2020
    
    """
       
    # set matlab path to SPM12 folder
    MatlabCommand.set_default_paths(pathSPM12)
    
    # parameters
    csf_max = 0.1 # c3 tissue class threshold.
    bone_max = 0.1 # c4 tissue class threshold.
    soft_max = 0.1 # c5 tissue class threshold.
    air_max = 0.1 # c6 tissue class threshold.
    
    # get path and file name
    path_split = os.path.split(filename)
    file = path_split[1]
    
    # make skullstrip folder
    path_skull = os.path.join(path_output, "skull")
    if not os.path.exists(path_skull):
        os.mkdir(path_skull)
        
    shutil.copyfile(filename, os.path.join(path_skull,file))
    os.chdir(path_skull)
    
    skull = NewSegment()
    skull.inputs.channel_files = os.path.join(path_skull,file)
    skull.inputs.channel_info = (0.001, 18, (True, True))
    skull.inputs.affine_regularization = "mni"
    skull.inputs.sampling_distance = 2
    skull.inputs.use_v8struct = True
    skull.inputs.warping_regularization = [0, 0.001, 0.5, 0.05, 0.2]
    skull.inputs.write_deformation_fields = [False, False]
    skull.inputs.mfile = True
    tissue1 = ((os.path.join(pathSPM12, "tpm/TPM.nii"), 1), 1, (True,False), (False, False))
    tissue2 = ((os.path.join(pathSPM12, "tpm/TPM.nii"), 2), 1, (True,False), (False, False))
    tissue3 = ((os.path.join(pathSPM12, "tpm/TPM.nii"), 3), 2, (True,False), (False, False))
    tissue4 = ((os.path.join(pathSPM12, "tpm/TPM.nii"), 4), 3, (True,False), (False, False))
    tissue5 = ((os.path.join(pathSPM12, "tpm/TPM.nii"), 5), 4, (True,False), (False, False))
    tissue6 = ((os.path.join(pathSPM12, "tpm/TPM.nii"), 6), 2, (True,False), (False, False))
    skull.inputs.tissues = [tissue1, tissue2, tissue3, tissue4, tissue5, tissue6]
    skull.run() 
    os.remove(os.path.join(path_skull,os.path.basename(filename)))
    
    # load tissue classes
    gm_img = nb.load(os.path.join(path_skull, "c1"+file))
    wm_img = nb.load(os.path.join(path_skull, "c2"+file))
    csf_img = nb.load(os.path.join(path_skull, "c3"+file))
    bone_img = nb.load(os.path.join(path_skull, "c4"+file))
    soft_img = nb.load(os.path.join(path_skull, "c5"+file))
    air_img = nb.load(os.path.join(path_skull, "c6"+file))
    
    gm_array = gm_img.get_fdata()
    wm_array = wm_img.get_fdata()
    csf_array = csf_img.get_fdata()
    bone_array = bone_img.get_fdata()
    soft_array = soft_img.get_fdata()
    air_array = air_img.get_fdata() 
    
    # generate skullstrip mask
    mask_array = gm_array + wm_array
    mask_array[mask_array > 0] = 1
    
    # get rid of pial noise
    mask_array[csf_array >= csf_max] = 0
    mask_array[bone_array >= bone_max] = 0
    mask_array[soft_array >= soft_max] = 0
    mask_array[air_array >= air_max] = 0
    
    # save skullstrip mask
    output = nb.Nifti1Image(mask_array, gm_img.affine, gm_img.header)
    nb.save(output,os.path.join(path_skull,"skullstrip_mask.nii"))
