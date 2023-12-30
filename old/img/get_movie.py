# -*- coding: utf-8 -*-

import os

import imageio as io
import nibabel as nb
import numpy as np


def get_movie(file_in, path_output, name_output, coord, axis=0, fps=10, transpose=True):
    """Get Movie.

    This function generates gif of a specific slice over a time series. It has
    the purpose to easily check preprocessing performance.

    Parameters
    ----------
    file_in : str
        File name of input file.
    path_output : str
        Path where output is saved.
    name_output : str
        Output file name without file extension.
    coord : int
        Selected slice.
    axis : int, optional
        Axis along which slice is selected. The default is 0.
    fps : float, optional
        Frames per second. The default is 10.
    transpose : bool, optional
        Rotate image. The default is True.

    Raises
    ------
    ValueError
        If `axis` number is invalid.

    Returns
    -------
    None.

    """

    # make sub-folder
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # load time series
    data_img = nb.load(file_in)
    data_array = data_img.get_fdata()

    # get matrix dimension
    x_size = data_img.header["dim"][1]
    y_size = data_img.header["dim"][2]
    z_size = data_img.header["dim"][3]
    nt = data_img.header["dim"][4]

    # get movie array
    if axis == 0:
        movie_array = np.zeros([y_size, z_size, nt])
        for i in range(nt):
            movie_array[:, :, i] = data_array[coord, :, :, i]
    elif axis == 1:
        movie_array = np.zeros([x_size, z_size, nt])
        for i in range(nt):
            movie_array[:, :, i] = data_array[:, coord, :, i]
    elif axis == 2:
        movie_array = np.zeros([x_size, y_size, nt])
        for i in range(nt):
            movie_array[:, :, i] = data_array[:, :, coord, i]
    else:
        raise ValueError("Choose a valid axis!")

    # scale movie array
    movie_array = movie_array / np.max(movie_array) * 255
    movie_array = movie_array.astype("uint8")

    images = []
    for i in range(nt):
        images.append(movie_array[:, :, i])

    if transpose is True:
        images = np.transpose(images, axes=(0, 2, 1))

    # generate video file
    io.mimwrite(os.path.join(path_output, name_output + ".gif"), images, fps=fps)
