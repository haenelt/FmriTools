# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb
from nipype.interfaces.fsl import ErodeImage, DilateImage
from scipy.ndimage.morphology import binary_fill_holes


def skullstrip_epi(file_in, roi_size=5, scale=0.75, nerode=2, ndilate=1,
                   savemask=False, cleanup=True):
    """Skullstrip EPI.

    Skullstrip input volume by defining an intensity threshold from the inner of 
    the brain volume. From a defined mid-point, a brain mask is grown inside the 
    brain. A binary filling holes algorithm is applied. To reduce remaining 
    skull within the brain mask, the mask is eroded and dilated several times.        

    Parameters
    ----------
    file_in : str
        Input file.
    roi_size : int, optional
        Size of cubic roi for image intensity threshold. The default is 5.
    scale : float, optional
        Scale image intensity threshold. The default is 0.75.
    nerode : int, optional
        Number of eroding iterations. The default is 2.
    ndilate : int, optional
        Number of dilating iterations. The default is 1.
    savemask : bool, optional
        Save mask time series. The default is False.
    cleanup : bool, optional
        Delete intermediate files after running. The default is True.

    Returns
    -------
    None.
    
    """
    
    # prepare path and filename
    path = os.path.split(file_in)[0]
    file = os.path.split(file_in)[1]

    # load data
    data_img = nb.load(os.path.join(path, file))
    data_array = data_img.get_data()
    
    # calculate mean intensity
    data_mean = data_array.mean()

    # get point within the brain
    inds = np.transpose(np.nonzero(data_array > data_mean))    
    x_mean = np.uint8(np.round((np.max(inds[:, 0])+np.min(inds[:, 0]))/2))
    y_mean = np.uint8(np.round((np.max(inds[:, 1])+np.min(inds[:, 1]))/2))
    z_mean = np.uint8(np.round((np.max(inds[:, 2])+np.min(inds[:, 2]))/2))
    
    # initialise mask
    mask_array = np.zeros_like(data_array)
    mask_temp_array = np.zeros_like(data_array)
    mask_array[x_mean, y_mean, z_mean] = 1
    mask_temp_array[x_mean, y_mean, z_mean] = 1
        
    # compute threshold
    roi = data_array[np.uint8(np.round(x_mean-roi_size/2)):np.uint8(np.round(x_mean+roi_size-1/2)),
                     np.uint8(np.round(y_mean-roi_size/2)):np.uint8(np.round(y_mean+roi_size-1/2)),
                     np.uint8(np.round(z_mean-roi_size/2)):np.uint8(np.round(z_mean+roi_size-1/2))]
    roi_mean = roi.mean()

    # grow mask
    while True:
            
        coords = np.transpose(np.nonzero(mask_array > 0))
    
        coords = coords[coords[:, 0] > 0, :]
        coords = coords[coords[:, 1] > 0, :]
        coords = coords[coords[:, 2] > 0, :]
        coords = coords[coords[:, 0] < np.size(data_array, 0)-1]
        coords = coords[coords[:, 1] < np.size(data_array, 1)-1]
        coords = coords[coords[:, 2] < np.size(data_array, 2)-1]
   
        # calculate neighbour coordinate
        mask_temp_array[coords[:, 0]-1, coords[:, 1], coords[:, 2]] = 1
        mask_temp_array[coords[:, 0], coords[:, 1]-1, coords[:, 2]] = 1
        mask_temp_array[coords[:, 0], coords[:, 1], coords[:, 2]-1] = 1
        mask_temp_array[coords[:, 0]+1, coords[:, 1], coords[:, 2]] = 1
        mask_temp_array[coords[:, 0], coords[:, 1]+1, coords[:, 2]] = 1
        mask_temp_array[coords[:, 0], coords[:, 1], coords[:, 2]+1] = 1
   
        # delete all old mask elements
        mask_temp_array = mask_temp_array - mask_array      
        mask_temp_array[data_array < scale*roi_mean] = 0
    
        # reinitialise mask_temp
        mask_array = mask_array + mask_temp_array

        # check break condition
        coords = np.transpose(np.nonzero(mask_temp_array == True))
        if len(coords) == 0:
            break

    # flood filling on brain mask
    mask_array = binary_fill_holes(mask_array, structure=np.ones((2, 2, 2)))

    # write mask (intermediate)
    newimg = nb.Nifti1Image(mask_array, data_img.affine, data_img.header)
    newimg.header['dim'][0] = 3
    nb.save(newimg, os.path.join(path, 'temp.nii'))

    # erode mask
    for i in range(nerode):
        erode = ErodeImage()
        erode.inputs.in_file = os.path.join(path, 'temp.nii')
        erode.inputs.output_type = 'NIFTI'
        erode.inputs.out_file = os.path.join(path, 'temp.nii')
        erode.run()

    # dilate mask
    for i in range(ndilate):
        dilate = DilateImage()
        dilate.inputs.in_file = os.path.join(path, 'temp.nii')
        dilate.inputs.operation = 'mean'
        dilate.inputs.output_type = 'NIFTI'
        dilate.inputs.out_file = os.path.join(path, 'temp.nii')
        dilate.run()
        
    # load final mask
    temp_img = nb.load(os.path.join(path, 'temp.nii'))
    mask_array = temp_img.get_data()
    
    # write masked image
    data_masked_array = data_array * mask_array
    output = nb.Nifti1Image(data_masked_array, data_img.affine, data_img.header)
    nb.save(output, os.path.join(path, 'p'+file))
    
    # write final output
    if savemask is True:
        newimg = nb.Nifti1Image(mask_array, data_img.affine, data_img.header)
        nb.save(newimg, os.path.join(path, 'mask_'+file))
    
    # clear output
    if cleanup is True:
        os.remove(os.path.join(path, 'temp.nii'))
