# -*- coding: utf-8 -*-

# python standard library inputs
import os
import shutil as sh

# external inputs
import numpy as np


def surface_flattening(fileREF, filePATCH, path_output, cleanup=True):
    """ Surface flattening
    
    Uses the FreeSurfer mris_flatten function to flatten a dense patch of 
    manually defined cortex. The manually cutted patch should have the following 
    file name <hemi>.<namePATCH>.patch.3d. Instead of smoothwm, we use white for 
    surface flattening.    

    Parameters
    ----------
    fileREF : str
        Reference surface file for flattening.
    filePATCH : str
        Path to be flattened saved as <hemi>.<namePATCH>.patch.3d.
    path_output : str
        Path where output is saved.
    cleanup : bool, optional
        Delete intermediate files. The default is True.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 01-11-2018             
    Last modified: 23-10-2020
    
    """
    
    # create temporary folder
    tmp = np.random.randint(0, 10, 5)
    tmp_string = ''.join(str(i) for i in tmp)
    path_temp = os.path.join(os.path.dirname(fileREF),"tmp_"+tmp_string)
    
    # make temporary folder
    if not os.path.exists(path_temp):
        os.mkdir(path_temp)
    else:
        raise FileExistsError("Temporary folder already exists!")

    # change to temporary folder
    cwd = os.getcwd()
    os.chdir(path_temp)
    
    # divide patch basename
    hemi = os.path.splitext(os.path.splitext(os.path.splitext(os.path.basename(filePATCH))[0])[0])[0]
    namePATCH = os.path.splitext(os.path.splitext(os.path.splitext(os.path.basename(filePATCH))[0])[0])[1]
    
    # copy reference file and path into temporary folder
    sh.copy2(fileREF,os.path.join(path_temp,hemi+".smoothwm"))
    sh.copy(filePATCH,os.path.basename(filePATCH))
    
    # surface flattening
    w = 0 # write out the surface every number of iterations.
    s = 20 # size of neighbourhood to be used in the optimization
    n = 7 # number of vertices at each distance to be used in the optimization
    os.system("mris_flatten" + \
              " -w " + str(w) + \
              " -distances " + str(s) + " " + str(n) + \
              " " + hemi + namePATCH + ".patch.3d" + \
              " " + hemi + namePATCH + ".patch.flat")
       
    # copy output
    sh.copy2(os.path.join(path_temp,hemi+namePATCH+".patch.flat"),
             os.path.join(path_output,hemi+namePATCH+".patch.flat"))
    sh.copy2(os.path.join(path_temp,hemi+namePATCH+".patch.flat.out"),
             os.path.join(path_output,hemi+namePATCH+".patch.flat.out"))
    
    # change to old folder
    os.chdir(cwd)
    
    # delete temporary files
    if cleanup:
        sh.rmtree(path_temp, ignore_errors=True)
