"""
Show low resolution of colour system

With this script, the low resolution of colour perception can be shown. An rgb image is read and 
converted to greyscale. Then, a grid is defined on the image and only the grid points are coloured
by a smoothed version of the initial colour image.

created by Daniel Haenelt
Date created: 29-07-2019             
Last modified: 29-07-2019  
"""
import os
import numpy as np
import imageio as io
from scipy.ndimage import gaussian_filter

input = "/home/daniel/Bilder/background.jpg"
path_output = "/home/daniel/Schreibtisch"

# parameters
filter_fwhm = 10 # gaussian smoothing of colour image
size_grid = 20 # grid size
thickness_grid = 5 # thickness of grid

""" do not edit below """

# load rgb input image
img = io.imread(input).astype(float)

# generate greyscale image
img_bw = np.sqrt(img[:,:,0]**2+img[:,:,1]**2+img[:,:,2]**2).astype(float)
img_bw = img_bw / np.max(img_bw) * 255

# compute grid coordinates
x = []
y = []
for i in range(thickness_grid):
    x = np.append(x,np.arange(i,np.shape(img_bw)[0],size_grid)).astype(int)
    y = np.append(y,np.arange(i,np.shape(img_bw)[1],size_grid)).astype(int)

# set grid points to zero
img_bw[x,:] = 0
img_bw[:,y] = 0

# gaussian smoothing of colour image
img_col = gaussian_filter(img, sigma=(filter_fwhm, filter_fwhm, 0), order=0)

# coloured grid image
img_grid = img_col.copy()
img_grid[:,:,0][img_bw != 0] = 0
img_grid[:,:,1][img_bw != 0] = 0
img_grid[:,:,2][img_bw != 0] = 0 

# final image with coloured grids
img_res = np.zeros_like(img)
img_res[:,:,0] = img_bw
img_res[:,:,1] = img_bw
img_res[:,:,2] = img_bw

img_res[x,:,0] = img_col[x,:,0]
img_res[x,:,1] = img_col[x,:,1]
img_res[x,:,2] = img_col[x,:,2]
img_res[:,y,0] = img_col[:,y,0]
img_res[:,y,1] = img_col[:,y,1]
img_res[:,y,2] = img_col[:,y,2]

io.imwrite(os.path.join(path_output,"img_bw.png"),img_bw.astype(int))
io.imwrite(os.path.join(path_output,"img_col.png"),img_col.astype(int))
io.imwrite(os.path.join(path_output,"img_grid.png"),img_grid.astype(int))
io.imwrite(os.path.join(path_output,"img_res.png"),img_res.astype(int))