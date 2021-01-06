# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys
import datetime
import subprocess
from shutil import copyfile

# external inputs
import numpy as np
from nibabel.freesurfer.io import read_geometry, write_geometry

# local inputs
from fmri_tools.io.get_filename import get_filename
from fmri_tools.surface.inflate_surf_mesh import inflate_surf_mesh


def _cart2pol(x, y, z):
    """
    Helper function for transformation cartesian to polar coordinates.
    """
    
    r = np.sqrt(x**2+y**2+z**2)
    phi = np.arctan2(y, x)
    theta = np.arccos(z/r)
    return r, phi, theta

    
def _pol2cart(r, phi, theta):
    """
    Helper function for transformation from polar to cartesian coordinates.
    """
    
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return x, y, z


def make_sphere(file_in, file_out, n_inflate=100, radius=None):
    """ Make sphere

    The scripts takes a generated FreeSurfer mesh and transformes it into
    a sphere with defined radius.    

    Parameters
    ----------
    file_in : str
        Filename of input surface.
    file_out : str
        Filename of output surface.
    n_inflate : int, optional
        Number of inflating iterations (if > 0). The default is 100.
    radius : float, optional
        Radius of final sphere in mm (if not None). The default is None.

    Returns
    -------
    None.

    None
    -------
    created by Daniel Haenelt
    Date created: 26-08-2020       
    Last modified: 25-10-2020    

    """
    
    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    
    # temporary file
    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = ''.join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    file_tmp = os.path.join(path_output, tmp_string)
    
    if os.path.exists(file_tmp):
        raise FileExistsError("Temporary file already exists!")

    # inflate surface mesh
    if n_inflate:
        inflate_surf_mesh(file_in, 
                          file_tmp, 
                          n_inflate)
    else:
        copyfile(file_in, file_tmp)
        
    
    # inflate surface
    try:
        subprocess.run(['mris_sphere', 
                        '-q', 
                        file_tmp, 
                        file_out], check = True)
    except subprocess.CalledProcessError:
        sys.exit("Sphere computation failed!")
    

    # change radius      
    if radius:
        vtx, fac = read_geometry(file_out)
        r, phi, theta = _cart2pol(vtx[:,0], vtx[:,1], vtx[:,2])
        r[:] = radius
        vtx[:,0], vtx[:,1], vtx[:,2] = _pol2cart(r, phi, theta)
        write_geometry(file_out, vtx, fac)
    
    # remove temporary file
    os.remove(file_tmp)
    