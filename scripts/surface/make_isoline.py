# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
from nibabel.freesurfer.io import read_geometry
from scipy.interpolate import griddata

# local inputs
from fmri_tools.mapping.map2grid import map2grid


"""
Make eccentricity isolines

This script computes isolines based on the phase map of a retinotopy 
measurement. The phase map is mapped onto the grid representation of a flattened 
surface patch. Only contour lines which are continuous are considered. The 
purpose of this script is to estimate cortical length D which is defined as 
isoline between dorsal and ventral V1/V2 borders. This can be used to define 
cortical frequency F = 1 / 2D similar to Parkes et al., 2005.
The script consists of the following steps:
    1. map phase onto grid
    2. smooth phase map
    3. get contour line
    4. check if line is continuous
    5. get cmap indices of line
    6. get white surface coordinates of cmap indices
    7. compute length

created by Daniel Haenelt
Date created: 15-02-2019
Last modified: 12-10-2020
"""

# input
input_white = "/data/pt_01880/V2STRIPES/p6/anatomy/dense/rh.white"
input_cmap = "/data/pt_01880/V2STRIPES/p6/anatomy/ortho/rh.occip4.patch.flat.cmap.nii"
input_mask = "/data/pt_01880/V2STRIPES/p6/anatomy/ortho/rh.occip4.patch.flat.mask.nii"
input_phase = "/data/pt_01880/V2STRIPES/p6/retinotopy/avg/surf/rh.ecc_phase_avg_def_layer5.mgh"
hemi = "rh"
name_patch = "occip4"
path_output = "/data/pt_01880/V2STRIPES/p6/anatomy/ortho"

# parameters
sigma = 20 # size of Gaussian filter
step_size = 10 # step size for contour lines
max_distance = 1 # maximal considered between vertex distance

# do not edit below

def gaussian_filter_edge(U, sigma):
    """
    This function calculates a Gaussian filter considering the voxels set to NaN 
    ouside the grid.
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

# load phase, cmap and mask
phase = map2grid(input_cmap, input_phase, 0, path_output, overwrite=False)
cmap = nb.load(input_cmap).get_fdata()
mask = nb.load(input_mask).get_fdata()
mask[mask==0] = np.NaN

# mask phase
phase = phase * mask

# apply gaussian filter to phase map
phase = gaussian_filter_edge(phase,sigma)
phase = phase * mask

# read white surface
coords, _ = read_geometry(input_white)

# get contour lines
D = np.zeros_like(phase)
D[D == 0] = np.NaN
for i in range(np.shape(phase)[0]):
    for j in range(np.shape(phase)[1]):
        
        # only consider pixels within the mask
        if np.isnan(phase[i,j]):
            continue
        else:
            line = plt.contour(phase,[phase[i,j]]).allsegs
                
        # discard not continuous line
        if np.size(line[0][0]) == np.size(line[0]):
            
            # get nearest pixels of contour line
            x = np.round(line[0][0][:,0]).astype(int)
            y = np.round(line[0][0][:,1]).astype(int)
            
            # edge consideration because of ceiling
            x[x >= np.shape(phase)[1]] = np.shape(phase)[1] - 1
            y[y >= np.shape(phase)[0]] = np.shape(phase)[0] - 1
            
            # get cmap indices of contour line
            ind = cmap[y,x].astype(int)
            
            # remove multiple indices
            ind_val, ind_ind = np.unique(ind,return_index=True)
            temp = np.zeros(len(ind))
            temp[:] = np.NaN
            temp[ind_ind] = ind_val
            ind = temp[~np.isnan(temp)].astype(int)
            
            # get white surface vertex coordinates of contour line
            x_white = coords[ind,0]
            y_white = coords[ind,1]
            z_white = coords[ind,2]
            
            # compute length between neighboring vertices
            delta_x = np.diff(x_white)
            delta_y = np.diff(y_white)
            delta_z = np.diff(z_white)
            
            # wrong cmap at the edges indices may lead to unrealistic lengths between neighbors 
            # which are discardeds
            if np.any( np.abs(np.concatenate((delta_x, delta_y, delta_z))) > max_distance ):
                continue
            
            # iso-eccentricity length from contour line
            D[i,j] = np.sum( np.sqrt(delta_x**2 + delta_y**2 + delta_z**2) )
        
        plt.close("all")

# get undefined pixel by linear interpolation
x, y = np.indices(D.shape)
D_interp = np.array(D)

D_interp[np.isnan(D_interp)] = griddata(
        (x[~np.isnan(D)], y[~np.isnan(D)]), # points we know
        D[~np.isnan(D)],                    # value we know
        (x[np.isnan(D)], y[np.isnan(D)]),
        "linear"
        )

# mask final distance map
D_interp = D_interp * mask

# initialize header
empty_header = nb.Nifti1Header()
empty_affine = np.eye(4)

# write distance map
output = nb.Nifti1Image(D_interp, empty_affine, empty_header)
nb.save(output, os.path.join(path_output,hemi+"."+name_patch+".iso_distance.nii"))

# Plot: Contours overlaid on grid
# load phase
phase_old = map2grid(input_cmap, input_phase, 0, path_output, overwrite=False)

# get contour lines
x1 = np.ceil(np.nanmin(phase))
x2 = np.floor(np.nanmax(phase))
lines = plt.contour(phase,np.arange(x1,x2,step_size)).allsegs

# remove not continuous lines
label = []
for i in range(len(lines)):
    if np.size(line[0][0]) == np.size(line[0]):
        label.append(np.concatenate(lines[i]))

# save plot
fig = plt.figure(frameon=False)
fig.set_size_inches(1,1)
ax = plt.Axes(fig,[0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
ax.imshow(phase_old, cmap="binary", aspect="auto")
for i in range(len(label)):
    ax.plot(label[i][:,0],label[i][:,1], linewidth=0.25, color="r")
fig.savefig(os.path.join(path_output,hemi+"."+name_patch+".isolines.png"), dpi=400)        
plt.close('all')
