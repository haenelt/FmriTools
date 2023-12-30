# -*- coding: utf-8 -*-
"""File name utilities."""

import os

__all__ = ["get_filename"]


def get_filename(file_in):
    """This function gets path, file name and file extension for an input filename. The
    loop checks for some given extension names. Otherwise it stops after the last found
    dot in the string.

    Parameters
    ----------
    file_in : str
        File name.

    Returns
    -------
    path : str
        Path of filename.
    name_file : str
        Basename of filename.
    ext_file : str
        File extension of filename.

    """
    # get path of input
    path = os.path.dirname(file_in)

    # get basename of input
    name_file = os.path.basename(file_in)

    # split basename and file extension
    ext_file = ""
    exit_loop = 0
    ext_key = [".nii", ".mgh", ".mgz"]
    while exit_loop == 0:
        name_file, ext_temp = os.path.splitext(name_file)
        ext_file = ext_temp + ext_file

        if not ext_temp:
            exit_loop = 1

        if ext_file in ext_key:
            exit_loop = 1

    return path, name_file, ext_file
