# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
from nibabel.freesurfer.io import read_morph_data, write_morph_data, read_geometry
from scipy.interpolate import griddata  


def morph2dense(source_sphere,target_sphere,input_morph,path_output):
    """ Morph to dense
    
    This function maps a morphological file from a source to a target surface.

    Parameters
    ----------
    source_sphere : str
        Source surface.
    target_sphere : str
        Target surface.
    input_morph : str
        Morphological input file.
    path_output : str
        Path where output is saved.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 13-07-2019
    Last modified: 12-10-2020

    """
    
    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # transform morphological data to dense surfaces
    pts_sphere_dense, _ = read_geometry(target_sphere)
    pts_sphere, _ = read_geometry(source_sphere)
    
    # get morphological data
    morph = read_morph_data(input_morph)
    
    # do the transformation
    method = "nearest"
    morph_dense = griddata(pts_sphere, morph, pts_sphere_dense, method)
        
    # write dense morphological data
    write_morph_data(os.path.join(path_output,os.path.basename(input_morph)), morph_dense)
    