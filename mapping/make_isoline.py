"""
Make eccentricity isolines

This script computes isolines based on the phase map of a retinotopy measurement. The phase map is
mapped onto the grid representation of a flattened surface patch. Only contour lines wich are 
continuous are considered. The script consists of the following steps:
    1. map phase onto grid
    2. smooth phase map
    3. get contour lines
    4. select continuous contour lines
    5. get cmap index of each vertex
    6. remove doubled indices
    7. save cmap indices as textfile

created by Daniel Haenelt
Date created: 15-02-2019
Last modified: 18-02-2019
"""
import os
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
from lib.mapping.map2grid import map2grid

# input
input_cmap = "/home/daniel/Schreibtisch/retinotopy/lh.occip4.patch.flat.cmap.nii"
input_mask = "/home/daniel/Schreibtisch/retinotopy/lh.occip4.patch.flat.mask.nii"
input_phase = "/home/daniel/Schreibtisch/retinotopy/lh.ecc_phase_avg_def_layer5.mgh"
path_output = "/home/daniel/Schreibtisch/bla"

# parameters
sigma = 20 # size of Gaussian filter
step_size = 10 # step size for contour lines

""" do not edit below """

def gaussian_filter_edge(U, sigma):
    """
    This function calculates a Gaussian filter considering the voxels set to NaN ouside the grid.
    """
    from scipy.ndimage.filters import gaussian_filter
    
    V = U.copy()
    V[U!=U] = 0
    VV = gaussian_filter(V,sigma)

    W = 0*U.copy()+1
    W[U!=U] = 0
    WW = gaussian_filter(W,sigma)

    return VV/WW

# make output folder
if not os.path.exists(path_output):
    os.mkdir(path_output)

# load phase map on grid
phase = map2grid(input_cmap, input_phase, 0, path_output, overwrite=False)
phase_old = phase.copy()

# load mask
mask = nb.load(input_mask).get_fdata()
mask[mask==0] = np.NaN

# mask phase mask
phase = phase * mask

# apply gaussian filter to phase map
phase = gaussian_filter_edge(phase,sigma)
phase = phase * mask

# get contour lines
x1 = np.ceil(np.nanmin(phase))
x2 = np.floor(np.nanmax(phase))
lines = plt.contour(phase,np.arange(x1,x2,step_size)).allsegs

# remove not continuous lines
label = []
for i in range(len(lines)):
    if not len(lines[i]) > 1:
        label.append(np.concatenate(lines[i]))

# transform contour lines to cmap indices
for i in range(len(label)):
    x = np.round(label[i][:,0]).astype(int)
    y = np.round(label[i][:,1]).astype(int)

    cmap = nb.load(input_cmap).get_fdata()
    data = cmap[y,x]        

    # remove multiple indices
    data_val, data_ind = np.unique(data,return_index=True)
    temp = np.zeros(len(data))
    temp[:] = np.NaN
    temp[data_ind] = data_val
    data = temp[~np.isnan(temp)]
    
    # save label
    np.savetxt(os.path.join(path_output,"isoline"+str(i)+".txt"),data)

# save plot
fig = plt.figure(frameon=False)
fig.set_size_inches(1,1)
ax = plt.Axes(fig,[0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(phase_old, cmap="binary", aspect="auto")
for i in range(len(label)):
    ax.plot(label[i][:,0],label[i][:,1], linewidth=0.25, color="r")
fig.savefig(os.path.join(path_output,"isolines.png"), dpi=400)        
plt.close('all')