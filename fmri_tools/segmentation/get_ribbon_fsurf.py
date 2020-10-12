# -*- coding: utf-8 -*-

# python standard library inputs
import os
import datetime

# external inputs
from nipype.interfaces.freesurfer import VolumeMask


def get_ribbon_fsurf(path, sub):
    """
    This function calculates the grey matter ribbon mask using a wrapped 
    freesurfer function. Old files are moved to the freesurfer trash folder 
    tagged with a date string as suffix.
    Inputs:
        *path: path to the freesurfer segmentation folder.
        *sub: name of the freesurfer segmentation folder.
        
    created by Daniel Haenelt
    Date created: 06-12-2018
    Last modified: 12-10-2020
    """

    # parameters
    hemi = ["lh","rh"]

    # label values for ribbon mask
    left_whitelabel = 2
    left_ribbonlabel = 3
    right_whitelabel = 41
    right_ribbonlabel = 42

    # output folder (freesurfer trash folder)
    path_trash = os.path.join(path,sub,"trash")
    
    # get date string for moved files
    date = datetime.datetime.now().strftime("%Y%m%d%H%M")

    # move old ribbon file into trash folder
    for i in range(len(hemi)):
        file_in = os.path.join(path,sub,"mri",hemi[i]+".ribbon.mgz")
        file_out = os.path.join(path_trash,hemi[i]+".ribbon_backup_"+date+".mgz")
        os.rename(file_in, file_out)

    file_in = os.path.join(path,sub,"mri","ribbon.mgz")
    file_out = os.path.join(path_trash,"ribbon_backup_"+date+".mgz")
    os.rename(file_in, file_out)

    # get ribbon file
    volmask = VolumeMask()
    volmask.inputs.left_whitelabel = left_whitelabel
    volmask.inputs.left_ribbonlabel = left_ribbonlabel
    volmask.inputs.right_whitelabel = right_whitelabel
    volmask.inputs.right_ribbonlabel = right_ribbonlabel
    volmask.inputs.lh_pial = os.path.join(path,sub,"surf","lh.pial")
    volmask.inputs.rh_pial = os.path.join(path,sub,"surf","rh.pial")
    volmask.inputs.lh_white = os.path.join(path,sub,"surf","lh.white")
    volmask.inputs.rh_white = os.path.join(path,sub,"surf","rh.white")
    volmask.inputs.subjects_dir = path
    volmask.inputs.subject_id = sub
    volmask.inputs.save_ribbon = True
    volmask.run()
