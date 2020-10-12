# -*- coding: utf-8 -*-

# python standard library inputs
import subprocess

# external inputs
import numpy as np
from numpy.linalg import inv


def vox2ras(file_in):
    """
    This function reads an input volume and computes the transformation between 
    voxel space and freesurfer vertex RAS (right-anterior-superior) coordinate 
    system from the header information. Transformation for both directions are 
    returned.
    Inputs:
        *file_in: filename of reference volume.
    Outputs:
        *vox2ras_tkr: vox2ras transformation matrix.
        *ras2voxs_tkr: ras2vox transformation matrix.
    
    created by Daniel Haenelt
    Date created: 18-12-2019
    Last modified: 12-10-2020
    """

    # get affine vox2ras-tkr and ras2vox-tkr transformation to reference volume
    transformation = subprocess.check_output(['mri_info', file_in, '--{}'.format("ras2vox-tkr")]).decode()
    
    # ignore if warning is stated in first line
    if transformation[:7] == "WARNING":
        i = 0
        while True:
            if transformation[i] == "\n":
                transformation = transformation[i+1:]
                break
            else:
                i += 1
    
    num_transformation = [[float(x) for x in line.split()] for line in transformation.split('\n') if len(line)>0]
    
    # get final transformation matriced as numpy array
    ras2vox_tkr = np.array(num_transformation)
    vox2ras_tkr = inv(np.array(num_transformation))
    
    return vox2ras_tkr, ras2vox_tkr
