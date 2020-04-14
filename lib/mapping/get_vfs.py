def get_vfs(input_sphere, input_white, input_patch, input_aparc, hemi, ecc_real, ecc_imag, pol_real, 
            pol_imag, path_output, fwhm_ecc = 4.0, fwhm_pol = 2.0, fwhm_vfs = 8.0, cleanup=True):
    """
    The purpose of the following function is to calculate the visual field sign map from retinotopy 
    data. The FREESURFER environment has to be set. The function was only tested with an annotation
    file not converted to the upsampled dense format. However, results seem to be correct and for
    the moment, nothing is changed. Different smoothing kernels can be set for eccentricity and
    polar angle. Eccentricity is much smoother and higher filters can be applied. The final visual 
    fieldsign map is saved as <hemi>.fieldsign in the output folder.
    Inputs:
        *input_sphere: input freesurfer sphere.
        *input_white: input white matter surface.
        *input_patch: flattened patch.
        *input_aparc: input annotation file <hemi>.aparc.annot
        *hemi: hemisphere.
        *ecc_real: real eccentricity values sampled on the surface.
        *ecc_imag: imaginary eccentricity values sampled on the surface.
        *pol_real: real polar angle values sampled on the surface.
        *pol_imag: imaginary polar angle values sampled on the surface.
        *path_output: path where output is saved.
        *fwhm_ecc: smoothing kernel for eccentricity input in mm.
        *fwhm_pol: smoothing kernel for polar angle input in mm.
        *fwhm_vfs: smoothing kernel for visual fieldsign calculation input in mm.
        *cleanup: delete intermediate files.

    created by Daniel Haenelt
    Date created: 11-12-2018             
    Last modified: 11-12-2018  
    """
    import os
    import numpy as np
    import shutil as sh

    # set freesurfer path environment
    os.environ["SUBJECTS_DIR"] = path_output

    # freesurfer subject
    tmp = np.random.randint(0, 10, 5)
    tmp_string = ''.join(str(i) for i in tmp)
    sub = "tmp_"+tmp_string

    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)

    # mimic freesurfer folder structure (with some additional folder for intermediate files)
    path_sub = os.path.join(path_output,sub)
    path_surf = os.path.join(path_sub,"surf")
    path_label = os.path.join(path_sub,"label")

    os.makedirs(path_sub)
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
