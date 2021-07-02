# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys

# external inputs
from nipype.interfaces.freesurfer import Curvature
    

def get_curvature(file_in, path_output, a=10):
    """Get curvature.
    
    This function calculates a curvature file for an input surface mesh using 
    freesurfer. The input file needs to have a prefix which indicates the 
    hemisphere of the surface mesh.    

    Parameters
    ----------
    file_in : str
        Filename of input surface.
    path_output : str
        Path where output is written.
    a : int, optional
        Number of smoothing iterations. The default is 10.

    Returns
    -------
    None.
    
    """
    
    # get hemi from filename
    hemi = os.path.splitext(os.path.basename(file_in))[0]
    if not hemi == "lh" and not hemi == "rh":
        sys.exit("Could not identify hemi from filename!")
    
    # calculate curvature file
    curv = Curvature()
    curv.inputs.in_file = file_in
    curv.inputs.save = True
    curv.inputs.averages = a  # for curv file
    curv.run()
    
    # rename mean curvature to curv
    os.rename(file_in+".H", os.path.join(path_output, hemi+".curv"))
    
    # delete gaussian curvature file and temporary input file
    os.remove(file_in+".K")
