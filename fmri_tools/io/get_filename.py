# -*- coding: utf-8 -*-

# python standard library inputs
import os


def get_filename(input):
    """
    This function gets path, file name and file extension for an input filename. 
    The loop checks for some given extension names. Otherwise it stops after the 
    last found dot in the string.
    Inputs:
        *input: filename.
    Outputs:
        *path: path of filename.
        *name_file: basename of filename.
        *ext_file: file extension of filename.
        
    created by Daniel Haenelt
    Date created: 09-12-2019   
    Last modified: 12-10-2020
    """

    # get path of input
    path = os.path.dirname(input)

    # get basename of input
    name_file = os.path.basename(input)
    
    # split basename and file extension
    ext_file = ""
    exit_loop = 0
    ext_key = [".nii",".mgh",".mgz"]
    while exit_loop == 0:
        name_file, ext_temp = os.path.splitext(name_file)
        ext_file = ext_temp + ext_file
        
        if not ext_temp:
            exit_loop = 1

        if ext_file in ext_key:
            exit_loop = 1

    return path, name_file, ext_file
