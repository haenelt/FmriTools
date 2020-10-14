# -*- coding: utf-8 -*-

# python standard library inputs
import os
from shutil import copyfile

# external inputs
import numpy as np
from sh import gunzip

# local inputs
from fmri_tools.io.get_filename import get_filename

    
def upsample_volume(file_in, file_out, dxyz=[0.4, 0.4, 0.4], rmode="Cu"):
    """ Upsample volume

    This function upsamples a nifti volume using the afni function 3dresample. 
    Before running the function, set the afni environment by calling AFNI in the 
    terminal. Output is an upsampled nifti volume.    

    Parameters
    ----------
    file_in : str
        Nifti input filename.
    file_out : str
        Nifti output filename.
    dxyz : list, optional
        Array of target resolution in single dimensions. The default is 
        [0.4, 0.4, 0.4].
    rmode : str, optional
        Interpolation methods (Linear, NN, Cu, Bk). The default is "Cu".

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 16-12-2019        
    Last modified: 12-10-2020
    
    """
    
    # get path and file extension of input file
    path_in, _, ext_in = get_filename(file_in)
    
    # make temporary copy of input file
    tmp = np.random.randint(0, 10, 5)
    tmp_string = ''.join(str(i) for i in tmp)
    file_tmp = os.path.join(path_in,"tmp_"+tmp_string+ext_in)
    copyfile(file_in, file_tmp)
    
    if os.path.splitext(file_tmp)[1] == ".gz":
        gunzip(file_tmp)
        file_tmp = os.path.splitext(file_tmp)[0]

    # upsample volume
    os.system("3dresample " + \
              "-dxyz " + str(dxyz[0]) + " " + str(dxyz[1]) + " " + str(dxyz[2]) + " " +\
              "-rmode " + str(rmode) + " " + \
              "-inset " + file_tmp + " " + \
              "-prefix " + file_out)
    
    # remove temporary copy
    os.remove(file_tmp)
    