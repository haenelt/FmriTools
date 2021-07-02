# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys
import subprocess

# local inputs
from ..io.get_filename import get_filename


def inflate_surf_mesh(file_in, file_out, n_iter):
    """Inflate surf mesh.

    The scripts takes a generated FreeSurfer surfaces and inflates it.    

    Parameters
    ----------
    file_in : str
        Filename of input surface.
    file_out : str
        Filename of output surface.
    n_iter : int
        Number of inflating iterations.

    Returns
    -------
    None.
    
    """
    
    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # inflate surface
    try:
        subprocess.run(['mris_inflate', 
                        '-n', str(n_iter), 
                        '-no-save-sulc', 
                        file_in, 
                        file_out], check=True)
    except subprocess.CalledProcessError:
        sys.exit("Surface inflation failed!")
