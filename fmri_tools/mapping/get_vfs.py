# -*- coding: utf-8 -*-

# python standard library inputs
import os
import datetime

# external inputs
import numpy as np
import shutil as sh


def get_vfs(input_sphere, input_white, input_patch, input_aparc, hemi, ecc_real, 
            ecc_imag, pol_real, pol_imag, path_output, fwhm_ecc=4.0, 
            fwhm_pol=2.0, fwhm_vfs=8.0, cleanup=True):
    """ Get VFS
    
    The purpose of the following function is to calculate the visual field sign 
    (vfs) map from retinotopy data. The FREESURFER environment has to be set. 
    The function was only tested with an annotation file not converted to the 
    upsampled dense format. However, results seem to be correct and for the 
    moment, nothing is changed. Different smoothing kernels can be set for 
    eccentricity and polar angle. Eccentricity is much smoother and higher 
    filters can be applied. The final visual fieldsign map is saved as 
    <hemi>.fieldsign in the output folder.

    Parameters
    ----------
    input_sphere : str
        Input freesurfer sphere.
    input_white : str
        Input white matter surface.
    input_patch : str
        Flattened patch.
    input_aparc : str
        Input annotation file <hemi>.aparc.annot
    hemi : str
        Hemisphere.
    ecc_real : str
        Real eccentricity values sampled on the surface.
    ecc_imag : str
        Imaginary eccentricity values sampled on the surface.
    pol_real : str
        Real polar angle values sampled on the surface.
    pol_imag : str
        Imaginary polar angle values sampled on the surface.
    path_output : str
        Path where output is saved.
    fwhm_ecc : float, optional
        Smoothing kernel for eccentricity input in mm. The default is 4.0.
    fwhm_pol : float, optional
        Smoothing kernel for polar angle input in mm. The default is 2.0.
    fwhm_vfs : float, optional
        Smoothing kernel for vfs calculation input in mm. The default is 8.0.
    cleanup : bool, optional
        Delete intermediate files.  The default is True.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 11-12-2018           
    Last modified: 25-10-2020 

    """

    # set freesurfer path environment
    os.environ["SUBJECTS_DIR"] = path_output

    # freesurfer subject
    tmp1 = np.random.randint(0, 10, 5)
    tmp1 = ''.join(str(i) for i in tmp1)
    tmp2 = datetime.datetime.now().strftime("%S%f")
    tmp_string = tmp1 + tmp2
    sub = "tmp_"+tmp_string

    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # mimic freesurfer folder structure (with some additional folder for intermediate files)
    path_sub = os.path.join(path_output,sub)
    path_surf = os.path.join(path_sub,"surf")
    path_label = os.path.join(path_sub,"label")

    if not os.path.exists(path_sub):
        os.makedirs(path_sub)
    else:
        raise FileExistsError("Temporary folder already exists!")

    os.makedirs(path_surf)
    os.makedirs(path_label)

    # copy surfaces to mimicked freesurfer folders
    sh.copyfile(input_white, os.path.join(path_surf,hemi+".white"))
    sh.copyfile(input_sphere, os.path.join(path_surf,hemi+".sphere"))
    sh.copyfile(input_patch, os.path.join(path_surf,hemi+".occip.patch.flat"))
    sh.copyfile(ecc_real, os.path.join(path_surf,hemi+".ecc_real.mgh"))
    sh.copyfile(ecc_imag, os.path.join(path_surf,hemi+".ecc_imag.mgh"))
    sh.copyfile(pol_real, os.path.join(path_surf,hemi+".pol_real.mgh"))
    sh.copyfile(pol_imag, os.path.join(path_surf,hemi+".pol_imag.mgh"))
    sh.copyfile(input_aparc, os.path.join(path_label,hemi+".aparc.annot"))

    # smooth surface
    input_data = ["ecc_real", "ecc_imag", "pol_real", "pol_imag"]
    input_fwhm = [fwhm_ecc, fwhm_ecc, fwhm_pol, fwhm_pol]
    for i in range(len(input_data)):
        os.system("mris_fwhm" + \
                  " --s " + sub + \
                  " --hemi " + hemi + \
                  " --smooth-only " + \
                  " --fwhm " + str(input_fwhm[i]) + \
                  " --i " + os.path.join(path_surf,hemi+"."+input_data[i]+".mgh") + \
                  " --o " + os.path.join(path_surf,hemi+"."+input_data[i]+"_smooth.mgh"))
    
    # file names of smoothed eccentricity and polar angle files
    ecc_real = os.path.join(path_surf,hemi+"."+input_data[0]+"_smooth.mgh")
    ecc_imag = os.path.join(path_surf,hemi+"."+input_data[1]+"_smooth.mgh")
    pol_real = os.path.join(path_surf,hemi+"."+input_data[2]+"_smooth.mgh")
    pol_imag = os.path.join(path_surf,hemi+"."+input_data[3]+"_smooth.mgh") 

    # create fiedsign map
    os.system("mri_fieldsign" + \
              " --s " + sub + \
              " --hemi " + hemi + \
              " --patch occip.patch.flat" + \
              " --new " + \
              " --eccen " + ecc_real + " " + ecc_imag + \
              " --polar "  + pol_real + " " + pol_imag + \
              " --fs " + os.path.join(path_output,hemi+".fieldsign.mgh") + \
              " --fwhm " + str(fwhm_vfs))

    # delete intermediate files
    if cleanup:
        sh.rmtree(path_sub, ignore_errors=True)
