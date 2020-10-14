# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import imageio as io


def get_gif(img_file, path_output, name_output, nsteps, duration):
    """ Get GIF

    Create gif movie from n input images with fading from one image to the next 
    image.    

    Parameters
    ----------
    img_file : list
        List of input images.
    path_output : str
        Path where output is saved.
    name_output : str
        Base name of output file.
    nsteps : int
        Number of generated transition images.
    duration : float
        Duration for each frame in seconds.

    Returns
    -------
    None.

    Notes
    -------
    created by Daniel Haenelt
    Date created: 17-11-2018
    Last modified: 12-10-2020

    """
    
    # append first list item to the end of the list
    img_file.append(img_file[0])

    images = []
    for i in range(len(img_file)-1):
    
        # input
        img1 = io.imread(img_file[i])
        img2 = io.imread(img_file[i+1])
    
        # generates a list with transitions from img1 to img2
        for j in range(nsteps):
            img = np.multiply(img1, (nsteps-j)/nsteps) + np.multiply(img2, j/nsteps)
            img = img.astype("uint8")
            images.append(img)

    # save output gif
    io.mimwrite(os.path.join(path_output,name_output+".gif"), images, duration=duration)
