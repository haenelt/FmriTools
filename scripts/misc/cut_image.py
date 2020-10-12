# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import imageio as io


"""
Cut image

In the following script, an image is cut to a target image size in px. In the 
final image, the outer borders are cut off.

created by Daniel Haenelt
Date created: 04-09-2019            
Last modified: 12-10-2020  
"""

input_image = "/home/raid2/haenelt/projects/paradigms/localiser/circle_in1.png"

# parameters
path_output = ""
name_output = "" # basename of output image
width = 768 # target width in px
height = 768 # target height in px

# do not edit below

# load image
img = io.imread(input_image)

# get original width and height in px
width_orig = np.shape(img)[1]
height_orig = np.shape(img)[0]

# get number of cut pixels in both direction
h_width = np.round(( width_orig - width ) / 2).astype(int)
h_height = np.round(( height_orig - width ) / 2).astype(int)

# get cut image
img = img[h_height:height_orig-h_height,h_width:width_orig-h_width]

# write output image
io.imwrite(os.path.join(path_output,name_output+".png"),img)
