# -*- coding: utf-8 -*-

# python standard library inputs
import os
import datetime
import subprocess
import shutil as sh
from shutil import copyfile

# external inputs
import numpy as np

# local inputs
from ..io.get_filename import get_filename


def intracortical_smoothing(file_surf, file_overlay, file_out, tan_size=0, 
                            rad_start=0, rad_size=1, tan_weights="gauss", 
                            cleanup=True):
    """Intracortical smoothing.
    
    This function applies the freesurfer function mris_smooth_intracortical
    (available since FreeSurfer 7) which enables simultaneous smoothing in 
    radial and tangential direction of the cortical sheet. [1]

    Parameters
    ----------
    file_surf : list[str]
        List of input surfaces which must be sorted from wm to pial.
    file_overlay : list[str]
        List of input overlays which must be sorted from wm to pial.
    file_out : str
        File name of smoothed output overlay.
    tan_size : float, optional
        Tangential extent of the smoothing kernel which is defined as the order
        of the vertex neighborhood (number of vertices; max = 6). The default 
        is 0.
    rad_start : int, optional
        Starting surface mesh of the intracortical smoothing kernel in radial
        direction (max = number of input surfaces). The default is 0 (white 
        surface).
    rad_size : int, optional
        Radial extent of the intracortical smoothing kernel (number of adjacent
        meshes; max = number of input surfaces). The default is 1 (no 
        smoothing).
    tan_weights : str, optional
        Weighting function for tangential smoothing. "gauss": tan_size = FWHM,
        "distance": tan_size = radius of the neighborhood around the central
        vertex of the smoothing kernel. The default is "gauss".
    cleanup : bool, optional
        Delete intermediate files. The default is True.

    Raises
    ------
    ValueError
        If `file_out` has an invalid file extension.
    FileExistsError
        If `path_temp` already exists.

    Returns
    -------
    None.

    References
    -------    
    .. [1] Blazejewska, AI, et al., Intracortical smoothing of small-voxel fMRI 
    data can provide increased detection power without spatial resolution 
    losses compared to conventional large-voxel fMRI data, NeuroImage 189, 
    601--614 (2019).

    """
   
    # check output file name
    path_output, name_output, ext_output = get_filename(file_out)
    if ext_output not in [".mgh", ".mgz"]:
        raise ValueError("Output file name is expected to have the file "
                         "extension mgh or mgz!")
    
    # make output folder
    create_folder = 0
    if not os.path.exists(path_output):
        create_folder = 1
        os.makedirs(path_output)
    
    # create temporary folder
    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = ''.join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    path_temp = os.path.join(path_output, "tmp_"+tmp_string)
    path_surf = os.path.join(path_temp, "surf")
    path_overlay = os.path.join(path_temp, "overlay")
    
    # make temporary folder
    if not os.path.exists(path_temp):
        os.mkdir(path_temp)
        os.mkdir(path_surf)
        os.mkdir(path_overlay)
    else:
        raise FileExistsError("Temporary folder already exists!")
    
    # copy input files to temporary folder
    for i, f in enumerate(file_surf):
        name_tmp = "surf_"+str(i).zfill(4)
        copyfile(f, os.path.join(path_surf, name_tmp))
    
    for i, f in enumerate(file_overlay):
        name_tmp = "overlay_"+str(i).zfill(4)+ext_output
        copyfile(f, os.path.join(path_overlay, name_tmp))

    # unix command
    command = [
        'mris_smooth_intracortical', 
        '--surf_dir', path_surf, 
        '--surf_name', 'surf_*',
        '--overlay_dir', path_overlay,
        '--overlay_name', 'overlay_*.mgh',
        '--output_dir', path_output,
        '--output_name', name_output+ext_output,
        '--tan-size', str(tan_size),
        '--rad-size', str(rad_size),
        '--rad-start', str(rad_start),
        '--tan-weights', str(tan_weights)]
    
    command = ' '.join(command)
    
    # run smoothing
    try:    
        subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError:
        sh.rmtree(path_temp, ignore_errors=True)
        if create_folder:
            sh.rmtree(path_output, ignore_errors=True)
        print("Intracortical smoothing failed!")
        
    # delete temporary files
    if cleanup:
        sh.rmtree(path_temp, ignore_errors=True)
