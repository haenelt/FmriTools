# -*- coding: utf-8 -*-

# python standard library inputs
import os
import shutil as sh

# external inputs
import nibabel as nb
from scipy.ndimage.morphology import binary_erosion
from nipype.interfaces.fsl import BET
from nipype.interfaces.freesurfer.preprocess import MRIConvert
from nighres.registration import apply_coordinate_mappings


def get_nuisance_mask(input, deformation, path_output, nerode_white=1, 
                      nerode_csf=1, segmentation=True, cleanup=True):
    """ Get nuisance mask
    
    This function calculates WM and CSF masks in space of the functional time 
    series. It uses SPM to compute WM and CSF probability maps. These maps are 
    masked with a skullstrip mask and transformed to native epi space.    

    Parameters
    ----------
    input : str
        Input anatomy (orig.mgz).
    deformation : str
        Coordinate mapping for ana to epi transformation.
    path_output : str
        Path where output is saved.
    nerode_white : int, optional
        Number of wm mask eroding steps. The default is 1.
    nerode_csf : int, optional
        Number of csf mask eroding steps. The default is 1.
    segmentation : bool, optional
        Do not calculate new masks to not rerun everything. The default is True.
    cleanup : bool, optional
        Delete intermediate files. The default is True.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 01-03-2019
    Last modified: 07-01-2021    

    """
    
    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # get filename without file extension of input file
    file = os.path.splitext(os.path.basename(input))[0]

    # convert to nifti format
    mc = MRIConvert()
    mc.inputs.in_file = input
    mc.inputs.out_file = os.path.join(path_output,file + ".nii")
    mc.inputs.out_type = "nii"
    mc.run()

    # bet skullstrip mask
    btr = BET()
    btr.inputs.in_file = os.path.join(path_output,file + ".nii")
    btr.inputs.frac = 0.5
    btr.inputs.mask = True
    btr.inputs.no_output = True
    btr.inputs.out_file = os.path.join(path_output,"bet")
    btr.inputs.output_type = "NIFTI"
    btr.run() 

    # segmentation
    if segmentation:
        os.system("matlab" + \
                  " -nodisplay -nodesktop -r " + \
                  "\"skullstrip_spm12(\'{0}\', \'{1}\'); exit;\"". \
                  format(os.path.join(path_output,file + ".nii"), 
                         path_output))

    # load tissue maps
    wm_array = nb.load(os.path.join(path_output,"skull","c2" + file + ".nii")).get_fdata()
    csf_array = nb.load(os.path.join(path_output,"skull","c3" + file + ".nii")).get_fdata()
    mask_array = nb.load(os.path.join(path_output,"bet_mask.nii")).get_fdata()

    # binarize
    wm_array[wm_array > 0] = 1
    csf_array[csf_array > 0] = 1

    # apply brain mask
    wm_array = wm_array * mask_array
    csf_array = csf_array * mask_array

    # erode wm
    wm_array = binary_erosion(
            wm_array, 
            structure=None, 
            iterations=nerode_white,
            mask=None, 
            output=None, 
            border_value=0, 
            origin=0, 
            brute_force=False,
            )

    # erode csf
    csf_array = binary_erosion(
            csf_array, 
            structure=None, 
            iterations=nerode_csf,
            mask=None, 
            output=None, 
            border_value=0, 
            origin=0, 
            brute_force=False,
            )

    # write wm and csf mask
    data_img = nb.load(input)
    wm_out = nb.Nifti1Image(wm_array, data_img.affine, data_img.header)
    nb.save(wm_out, os.path.join(path_output,"wm_mask_orig.nii"))
    csf_out = nb.Nifti1Image(csf_array, data_img.affine, data_img.header)
    nb.save(csf_out, os.path.join(path_output,"csf_mask_orig.nii"))

    # apply deformation to mask
    apply_coordinate_mappings(os.path.join(path_output,"wm_mask_orig.nii"), # input 
                              deformation, # cmap
                              interpolation = "nearest", # nearest or linear
                              padding = "zero", # closest, zero or max
                              save_data = True, # save output data to file (boolean)
                              overwrite = True, # overwrite existing results (boolean)
                              output_dir = path_output, # output directory
                              file_name = "wm_mask" # base name with file extension for output
                              )

    apply_coordinate_mappings(os.path.join(path_output,"csf_mask_orig.nii"), # input 
                              deformation, # cmap
                              interpolation = "nearest", # nearest or linear
                              padding = "zero", # closest, zero or max
                              save_data = True, # save output data to file (boolean)
                              overwrite = True, # overwrite existing results (boolean)
                              output_dir = path_output, # output directory
                              file_name = "csf_mask" # base name with file extension for output
                              )
    
    # rename transformed masks
    os.rename(os.path.join(path_output,"wm_mask_def-img.nii.gz"),
              os.path.join(path_output,"wm_mask.nii.gz"))
    os.rename(os.path.join(path_output,"csf_mask_def-img.nii.gz"),
              os.path.join(path_output,"csf_mask.nii.gz"))

    # cleanup
    if cleanup:
        os.remove(os.path.join(path_output,"bet_mask.nii"))
        os.remove(os.path.join(path_output,"csf_mask_orig.nii"))
        os.remove(os.path.join(path_output,"wm_mask_orig.nii"))
        os.remove(os.path.join(path_output,"orig.nii"))
        sh.rmtree(os.path.join(path_output,"skull"), ignore_errors=True)
