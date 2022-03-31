# -*- coding: utf-8 -*-

# python standard library inputs
import os
import sys
import datetime
import shutil as sh

# external inputs
import numpy as np
from nibabel.freesurfer.io import read_geometry
from nipype.interfaces.freesurfer import SampleToSurface

# local inputs
from ..io.get_filename import get_filename
from ..io.surf import write_mgh, read_mgh
from ..io.mgh2nii import mgh2nii


def map2surface(input_surf, input_vol, write_output=False, path_output="",
                interp_method="nearest", input_surf_target=None, input_ind=None,
                cleanup=True):
    """Map to surface.

    This function samples data from the input volume to the input surface and 
    optionally maps those values to a target surface if an index file is given.    

    Parameters
    ----------
    input_surf : str
        Surface mesh onto which volume data is sampled.
    input_vol : str
        Volume from which data is sampled.
    write_output : bool, optional
        Write sampled data as MGH file. The default is False.
    path_output : str, optional
        Path where to save output. The default is "".
    interp_method : str, optional
        Interpolation method (nearest or trilinear). The default is "nearest".
    input_surf_target : str, optional
        Target surface (only necessary if index file is given). The default is 
        None.
    input_ind : str, optional
        Textfile with mapping of vertex indices to target space. The default is 
        None.
    cleanup : bool, optional
        Remove intermediate files. The default is True.

    Raises
    ------
    FileExistsError
        If the temporary folder already exists.

    Returns
    -------
    arr_sampled : ndarray
        Image array.
    affine_sampled : ndarray
        Affine transformation matrix.
    header_sampled : MGHHeader
        Image header.

    """

    # clean everything if no output is written
    if not write_output:
        path_output, _, _ = get_filename(input_vol)
        #cleanup = True

    # set freesurfer path environment
    os.environ["SUBJECTS_DIR"] = path_output

    # freesurfer subject
    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = ''.join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    sub = "tmp_" + tmp_string

    # make output folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # mimic freesurfer folder structure (with some additional folder for 
    # intermediate files)
    path_sub = os.path.join(path_output, sub)
    path_mri = os.path.join(path_sub, "mri")
    path_surf = os.path.join(path_sub, "surf")

    if not os.path.exists(path_sub):
        os.makedirs(path_sub)
    else:
        raise FileExistsError("Temporary folder already exists!")

    os.makedirs(path_mri)
    os.makedirs(path_surf)

    # get filenames
    _, name_vol, ext_vol = get_filename(input_vol)
    _, hemi, name_surf = get_filename(input_surf)
    name_surf = name_surf.replace(".", "")

    # check filename
    if not hemi == "lh" and not hemi == "rh":
        sys.exit("Could not identify hemi from filename!")

    # copy input volume as orig.mgz to mimic freesurfer folder
    if ext_vol != ".mgz":
        mgh2nii(input_vol, path_mri, "mgz")
        os.rename(os.path.join(path_mri, name_vol + ".mgz"),
                  os.path.join(path_mri, "orig.mgz"))
    else:
        sh.copyfile(input_vol,
                    os.path.join(path_mri, "orig.mgz"))

    # copy input surface to mimic freesurfer folder
    sh.copyfile(input_surf, os.path.join(path_surf, hemi + ".source"))

    # filename of sampled data
    file_sampled = os.path.join(path_surf, hemi + "." + "sampled.mgh")

    # mri_vol2surf
    sampler = SampleToSurface()
    sampler.inputs.subject_id = sub
    sampler.inputs.reg_header = True
    sampler.inputs.hemi = hemi
    sampler.inputs.source_file = input_vol
    sampler.inputs.surface = "source"
    sampler.inputs.sampling_method = "point"
    sampler.inputs.sampling_range = 0
    sampler.inputs.sampling_units = "mm"
    sampler.inputs.interp_method = interp_method  # nearest or trilinear
    sampler.inputs.out_type = "mgh"
    sampler.inputs.out_file = file_sampled
    sampler.run()

    # load data
    arr_sampled, affine_sampled, header_sampled = read_mgh(file_sampled)

    # map on separate mesh
    if input_ind:
        # load data
        ind_target = np.loadtxt(input_ind, dtype=int)
        vtx_target, _ = read_geometry(input_surf_target)

        # read sampled morph data
        arr_tmp = arr_sampled.copy()

        # update header
        header_sampled["dims"][0] = len(vtx_target)
        header_sampled["Mdc"] = np.eye(3)

        # sample array in target space
        arr_sampled = np.zeros(len(vtx_target))
        arr_sampled[ind_target] = arr_tmp

    if write_output:
        file_out = os.path.join(path_output,
                                hemi + "." + name_vol + "_" + name_surf)

        if input_ind:
            file_out += "_trans.mgh"
        else:
            file_out += ".mgh"

        write_mgh(file_out,
                  arr_sampled,
                  affine_sampled,
                  header_sampled)

        # delete intermediate files
    if cleanup:
        sh.rmtree(path_sub, ignore_errors=True)
        if not len(os.listdir(path_output)):
            sh.rmtree(path_output)

    return arr_sampled, affine_sampled, header_sampled
