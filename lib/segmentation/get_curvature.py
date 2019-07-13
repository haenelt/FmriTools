def get_curvature(path,sub):
    """
    This function calculates a curvature file using freesurfer. Old files are moved to the  
    freesurfer trash folder tagged with a date string as suffix.
    Inputs:
        *path: path to the freesurfer segmentation folder.
        *sub: name of the freesurfer segmentation folder.
        
    created by Daniel Haenelt
    Date created: 13-07-2019
    Last modified: 13-07-2019
    """
    import os
    import datetime
    import numpy as np
    from nipype.interfaces.freesurfer import Curvature
   
    # parameters 
    hemi = ["lh","rh"] # hemisphere prefix
    file = "white" # surface from which curvature is calculated
    averages = 10 # number of smoothing iterations
       
    # output folder (freesurfer trash folder)
    path_trash = os.path.join(path,sub,"trash")
    
    # get date string for moved files
    date = datetime.datetime.now().strftime("%Y%m%d%H%M")

    # rename input file to prevent any overwriting
    tmp = np.random.randint(0, 10, 5)
    tmp_string = ''.join(str(i) for i in tmp)

    # move old surfaces and morphological files to trash folder
    for i in range(len(hemi)):
        file_in = os.path.join(path,sub,"surf",hemi[i]+".curv")
        file_out = os.path.join(path_trash,hemi[i]+".curv_backup_"+date)
        os.rename(file_in, file_out)

    # calculate curvature file
    for i in range(len(hemi)):
        
        os.rename(os.path.join(path,"surf",hemi+"."+file),
                  os.path.join(path,"surf",hemi+"."+tmp_string))
        
        curv = Curvature()
        curv.inputs.in_file = os.path.join(path,"surf",hemi+"."+tmp_string)
        curv.inputs.save = True
        curv.inputs.averages = averages # for curv file
        curv.run()
        
        # rename mean curvature to curv
        os.rename(os.path.join(path,"surf",hemi+"."+tmp_string+".H"),
                  os.path.join(path,"surf",hemi+".curv"))
        
        # delete gaussian curvature file and temporary input file
        os.remove(os.path.join(path,"surf",hemi+"."+tmp_string))
        os.remove(os.path.join(path,"surf",hemi+"."+tmp_string+".K"))