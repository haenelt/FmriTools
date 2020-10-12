# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb
import imageio as io
    

def get_movie(input, path_output, name_output, coord, axis=0, fps = 10, 
              transpose=True):
    """
    This function generates gif of a specific slice over a time series. It has 
    the purpose to easily check preprocessing performance.
    Inputs:
        *input: input file.
        *path_output: path where output is saved.
        *name_output: output file name without file extension.
        *coord: selected slice.
        *axis: axis along which slice is selected.
        *fps: frames per second.
        *transpose: rotate image.
        
    created by Daniel Haenelt
    Date created: 04-02-2019
    Last modified: 12-10-2020
    """
    
    # make subfolders
    if not os.path.exists(path_output):
        os.makedirs(path_output)
    
    # load time series
    data_img = nb.load(input)
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
            movie_array[:,:,i] = data_array[coord,:,:,i]
    elif axis == 1:
        movie_array = np.zeros([x_size, z_size, nt])
        for i in range(nt):
            movie_array[:,:,i] = data_array[:,coord,:,i]
    elif axis == 2:
        movie_array = np.zeros([x_size,y_size,nt])
        for i in range(nt):
            movie_array[:,:,i] = data_array[:,:,coord,i]
    else:
        print("Choose a valid axis!")
    
    # scale movie array
    movie_array = movie_array / np.max(movie_array) * 255
    movie_array = movie_array.astype('uint8')
    
    images = []
    for i in range(nt):
        images.append(movie_array[:,:,i])
    
    if transpose is True:
        images = np.transpose(images, axes = (0,2,1))
    
    # generate video file
    io.mimwrite(os.path.join(path_output,name_output+".gif"), images, fps = fps)
