# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys
import subprocess

# local inputs
from ..io.get_filename import get_filename


def smooth_surface(file_in, file_out, n_iter):
    """Smooth surface.
    
    This function smoothes a surface mesh using freesurfer.    
    
    Parameters
    ----------
    file_in : str
        Filename of input surface.
    file_out : str
        Filename of output surface.
    n_iter : int
        Number of smoothing iterations.

    Returns
    -------
    None.

    """
       
    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # smooth surface
    try:
        subprocess.run(['mris_smooth', 
                        '-n', str(n_iter), 
                        '-nw', 
                        file_in, 
                        file_out], check=True)
    except subprocess.CalledProcessError:
        sys.exit("Surface smoothing failed!")
