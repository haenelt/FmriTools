# -*- coding: utf-8 -*-

# python standard library inputs
import os
import shutil as sh

# external inputs
import numpy as np
import nibabel as nb
from nipype.interfaces.fsl import ConvertWarp, Merge
from nipype.interfaces.fsl.maths import MeanImage
from nipype.interfaces.fsl.preprocess import ApplyWarp

# local inputs
from ..io.get_filename import get_filename
from ..cmap.generate_coordinate_mapping import generate_coordinate_mapping


def gnl_correction(file_in, file_bash, file_coeff, python3_env, python2_env,
                   path_output, cleanup=True):
    """GNL correction.

    The purpose of the following function is to correct for gradient 
    nonlinearities. A corrected file is written using spline interpolation. The 
    function needs FSL to be included in the search path.    

    Parameters
    ----------
    file_in : str
        Filename of input image.
    file_bash : str
        Filename of bash script which calls the gradient unwarping toolbox.
    file_coeff : str
        Filename of siemens coefficient file.
    python3_env : str
        Name of python3 virtual environment.
    python2_env : str
        Name of python2 virtual environment.
    path_output : str
        Path where output is written.
    cleanup : bool, optional
        Delete intermediate files. The default is True.

    Returns
    -------
    None.

    """

    # get fileparts
    path, name, ext = get_filename(file_in)

    # make subfolders
    path_grad = os.path.join(path_output, "grad")
    if not os.path.exists(path_grad):
        os.makedirs(path_grad)

    # parse arguments
    file_output = os.path.join(path_output, name + "_gnlcorr" + ext)
    file_warp = os.path.join(path_grad, "warp.nii.gz")
    file_jacobian = os.path.join(path_grad, "warp_jacobian.nii.gz")
    
    # run gradient unwarp
    os.system("bash " + file_bash +
              " " + python3_env +
              " " + python2_env +
              " " + path_grad +
              " " + file_in +
              " trilinear.nii.gz" +
              " " + file_coeff)

    # now create an appropriate warpfield output (relative convention)
    convertwarp = ConvertWarp()
    convertwarp.inputs.reference = os.path.join(path_grad, "trilinear.nii.gz")
    convertwarp.inputs.warp1 = os.path.join(path_grad, "fullWarp_abs.nii.gz")
    convertwarp.inputs.abswarp = True
    convertwarp.inputs.out_relwarp = True
    convertwarp.inputs.out_file = file_warp
    convertwarp.inputs.args = "--jacobian=" + file_jacobian
    convertwarp.run()

    # convertwarp's jacobian output has 8 frames, each combination of one-sided
    # differences, so average them
    calcmean = MeanImage()
    calcmean.inputs.in_file = file_jacobian
    calcmean.inputs.dimension = "T"
    calcmean.inputs.out_file = file_jacobian
    calcmean.run()

    # apply warp to first volume
    applywarp = ApplyWarp()
    applywarp.inputs.in_file = file_in
    applywarp.inputs.ref_file = file_in
    applywarp.inputs.relwarp = True
    applywarp.inputs.field_file = file_warp
    applywarp.inputs.output_type = "NIFTI"
    applywarp.inputs.out_file = file_output
    applywarp.inputs.interp = "spline"
    applywarp.run()

    # normalise warped output image to initial intensity range
    data_img = nb.load(file_in)
    data_array = data_img.get_fdata()
    max_data = np.max(data_array)
    min_data = np.min(data_array)

    data_img = nb.load(file_output)
    data_array = data_img.get_fdata()
    data_array[data_array < min_data] = 0
    data_array[data_array > max_data] = max_data

    output = nb.Nifti1Image(data_array, data_img.affine, data_img.header)
    nb.save(output, file_output)

    # calculate gradient deviations
    os.system("calc_grad_perc_dev" +
              " --fullwarp=" + file_warp +
              " -o " + os.path.join(path_grad, "grad_dev"))

    # merge directions
    merger = Merge()
    merger.inputs.in_files = [os.path.join(path_grad, "grad_dev_x.nii.gz"),
                              os.path.join(path_grad, "grad_dev_y.nii.gz"),
                              os.path.join(path_grad, "grad_dev_z.nii.gz")]
    merger.inputs.dimension = 't'
    merger.inputs.merged_file = os.path.join(path_grad, "grad_dev.nii.gz")
    merger.run()

    # convert from % deviation to absolute
    data_img = nb.load(os.path.join(path_grad, "grad_dev.nii.gz"))
    data_array = data_img.get_fdata()
    data_array = data_array / 100

    output = nb.Nifti1Image(data_array, data_img.affine, data_img.header)
    nb.save(output, os.path.join(path_grad, "grad_dev.nii.gz"))

    # warp coordinate mapping
    generate_coordinate_mapping(file_in, 0, path_grad, suffix="gnl", time=False,
                                write_output=True)

    applywarp = ApplyWarp()
    applywarp.inputs.in_file = os.path.join(path_grad, "cmap_gnl.nii")
    applywarp.inputs.ref_file = file_in
    applywarp.inputs.relwarp = True
    applywarp.inputs.field_file = file_warp
    applywarp.inputs.out_file = os.path.join(path_grad, "cmap_gnl.nii")
    applywarp.inputs.interp = "trilinear"
    applywarp.inputs.output_type = "NIFTI"
    applywarp.run()

    # clean intermediate files
    if cleanup:
        sh.rmtree(path_grad, ignore_errors=True)
