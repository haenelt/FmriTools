# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb


def estimate_pv(input_target, input_border, path_output, name_output):
    """ Estimate PV

    This function estimates the partial volume contribution in each image voxel 
    of a target image from an upsampled binary image depicting the GM/WM or 
    GM/CSF border. Partial voluming is estimated by downsampling the binary 
    image using a moving-average like algorithm and calculating the ratio of 
    both binary elements within each target voxel.    

    Parameters
    ----------
    input_target : str
        Target space for which partial voluming is estimated.
    input_border : str
        Upsampled binary border depicting the high-resolution tissue border.
    path_output : str
        Path where output is saved.
    name_output : str
        Basename of output image.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 27-04-2019
    Last modified: 12-10-2020
    
    """
    
    # load data
    target = nb.load(input_target)
    border = nb.load(input_border)
    border_array = border.get_fdata()

    # downsampling parameters
    matrix_up = border.header["dim"][1:4]
    matrix_down = target.header["dim"][1:4]
    magn_factor = matrix_up/matrix_down

    # sampling steps
    p = np.arange(0,matrix_down[0]*magn_factor[0],magn_factor[0])
    q = np.arange(0,matrix_down[1]*magn_factor[1],magn_factor[1])
    r = np.arange(0,matrix_down[2]*magn_factor[2],magn_factor[2])

    # initialise matrices
    M = np.zeros(target.header["dim"][1:4]);
    start_weight = np.zeros(3)
    end_weight = np.zeros(3)

    # downsampling using a moving-average like algorithm
    for i in range(len(p)):
           
        # weightings (start_weight)
        if np.mod(p[i],1) != 0:
            start_weight[0] = 1 - np.mod(p[i],1)
        else:
            start_weight[0] = 1
        
        # weightings (end_weight)
        if i == len(p)-1:
            end_weight[0] = 1
        else:
            if np.mod(p[i+1],1) != 0:
                end_weight[0] = np.mod(p[i+1],1)
            else:
                end_weight[0] = 1
        
        for j in range(len(q)):
            
            # weightings (start_weight)
            if np.mod(q[j],1) != 0:
                start_weight[1] = 1 - np.mod(q[j],1)
            else:
                start_weight[1] = 1
        
            # weightings (end_weight)
            if j == len(q)-1:
                end_weight[1] = 1
            else:
                if np.mod(q[j+1],1) != 0:
                    end_weight[1] = np.mod(q[j+1],1)
                else:
                    end_weight[1] = 1
            
            for k in range(len(r)):
                
                # weightings (start_weight)
                if np.mod(r[k],1) != 0:
                    start_weight[2] = 1 - np.mod(r[k],1)
                else:
                    start_weight[2] = 1
        
                # weightings (end_weight)
                if k == len(r)-1:
                    end_weight[2] = 1
                else:
                    if np.mod(r[k+1],1) != 0:
                        end_weight[2] = np.mod(r[k+1],1)
                    else:
                        end_weight[2] = 1
                
                # p-coordinates
                if i == len(p)-1:
                    x = np.arange(np.round(p[-1]-0.5),matrix_up[0],1).astype(int)
                else:
                    x = np.arange(np.round(p[i]-0.5),np.round(p[i+1]-0.5),1).astype(int)
                
                # q-coordinates
                if j == len(q)-1:
                    y = np.arange(np.round(q[-1]-0.5),matrix_up[1],1).astype(int)
                else:
                    y = np.arange(np.round(q[j]-0.5),np.round(q[j+1]-0.5),1).astype(int)
                
                # r-coordinates
                if k == len(r)-1:
                    z = np.arange(np.round(r[-1]-0.5),matrix_up[2],1).astype(int)
                else:
                    z = np.arange(np.round(r[k]-0.5),np.round(r[k+1]-0.5),1).astype(int)
                
                # define submatrix
                X = border_array[x,:,:]
                X = X[:,y,:]
                X = X[:,:,z]
                
                # apply weightings
                X[0,:,:] = start_weight[0] * X[0,:,:];
                X[:,0,:] = start_weight[1] * X[:,0,:];
                X[:,:,0] = start_weight[2] * X[:,:,0];
                X[-1,:,:] = end_weight[0] * X[-1,:,:];
                X[:,-1,:] = end_weight[1] * X[:,-1,:];
                X[:,:,-1] = end_weight[2] * X[:,:,-1];                   
                
                X = np.reshape(X, np.size(X))
                X = X[~np.isnan(X)]
                
                # voxel value in target space
                M[i,j,k] = np.sum(X) / len(X)
    
    # save data
    output = nb.Nifti1Image(M, target.affine, target.header)
    nb.save(output, os.path.join(path_output,name_output+"_pve.nii"))
    