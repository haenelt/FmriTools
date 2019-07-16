def smooth_surface(path,sub,surf,niter):
    """
    This function smoothes a surface file using freesurfer. Old files are moved to the freesurfer 
    trash folder tagged with a date string as suffix.
    Inputs:
        *path: path to the freesurfer segmentation folder.
        *sub: name of the freesurfer segmentation folder.
        *surf: basename of surface file.
        *niter: number of smoothing iterations.
        
    created by Daniel Haenelt
    Date created: 13-07-2019
    Last modified: 13-07-2019
    """
    import os
    import datetime
    from nipype.interfaces.freesurfer import SmoothTessellation
   
    # parameters 
    hemi = ["lh","rh"] # hemisphere prefix
       
    # output folder (freesurfer trash folder)
    path_trash = os.path.join(path,sub,"trash")
    
    # get date string for moved files
    date = datetime.datetime.now().strftime("%Y%m%d%H%M")

    # move old surfaces and morphological files to trash folder
    for i in range(len(hemi)):
        file_in = os.path.join(path,sub,"surf",hemi[i]+"."+surf)
        file_out = os.path.join(path_trash,hemi[i]+"."+surf+"_backup_"+date)
        os.rename(file_in, file_out)

    # calculate curvature file
    for i in range(len(hemi)):
        smooth = SmoothTessellation()
        smooth.inputs.in_file = os.path.join(path_trash,hemi[i]+"."+surf+"_backup_"+date)
        smooth.inputs.out_file = os.path.join(path,sub,"surf",hemi[i]+"."+surf)
        smooth.inputs.smoothing_iterations = niter
        smooth.inputs.disable_estimates = True
        smooth.run()