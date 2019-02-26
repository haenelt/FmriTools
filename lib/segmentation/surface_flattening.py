def surface_flattening(path,path_dense,sub,hemi,namePATCH):
    """
    Uses the FreeSurfer mris_flatten function to flatten a dense patch of manually defined cortex. 
    The manually cutted patch should have the following file name <hemi>.<namePATCH>.patch.3d and 
    should be saved in the dense folder. 
    Inputs:
        *path: path to SUBJECTS_DIR of the freesurfer segmentation.
        *path_dense: path to the dense <hemi>.smoothwm and to <hemi>.<namePATCH>.patch.3d.
        *sub: freesurfer subject name.
        *hemi: hemisphere.
        *namePATCH: name of the cutted patch.

    created by Daniel Haenelt
    Date created: 01-11-2018             
    Last modified: 17-12-2018
    """
    import os
    import shutil
  
    # change directory
    os.chdir(path_dense)

    # copy dense smoothwm into freesurfer surf folder and temporary rename original smoothwm
    smoothwm_orig = os.path.join(path,sub,"surf",hemi+".smoothwm")
    smoothwm_dense = os.path.join(path_dense,hemi+".smoothwm")
    
    shutil.copy2(smoothwm_orig,smoothwm_orig+"_temp")
    os.remove(smoothwm_orig)
    shutil.copy2(smoothwm_dense,smoothwm_orig)
    
    # surface flattening
    w = 0 # write out the surface every number of iterations.
    s = 20 # size of neighbourhood to be used in the optimization
    n = 7 # number of vertices at each distance to be used in the optimization
    os.system("mris_flatten" + \
              " -w " + str(w) + \
              " -distances " + str(s) + " " + str(n) + \
              " " + hemi + "." + namePATCH + ".patch.3d" + \
              " " + hemi + "." + namePATCH + ".patch.flat")
    
    # remove dense smoothwm from freesurfer surf folder
    os.remove(smoothwm_orig)
    os.rename(smoothwm_orig+"_temp",smoothwm_orig)
