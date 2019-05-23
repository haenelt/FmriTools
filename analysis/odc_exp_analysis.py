"""
Analysis of experimental ocular dominance columns

The purpose of the following script is to analysis ocular dominance columns and estimate their 
column width and column spacing.

To do so, the power spectrum and the normalized autocorrelation (NAC) map are computed and evaluated 
along different projection lines. Especially, principal axes are first estimated from thresholded 
Fourier spectograms and width and spacing. 

In total, the following quantities are computed:

(1) column width is estimated by computing the FWHM of the NAC central peak.
(2) column spacing is estimated by computing the distance from the NAC central peak to its first 
neighbor.
(3) take amplitude of first neighbor peak as repetition marker?
(3) column spacing is further estimated by computing the spatial frequency of the maximum Fourier
component.
(4) take the maximum power (relative to central peak) as repetition marker?

created by Daniel Haenelt
Date created: 15-04-2019
Last modified: 17-04-2019
"""
import os
import numpy as np
import nibabel as nb
from lib.analysis import get_pca
from lib.analysis import analyze_fft
from lib.analysis import analyze_acorr

# parameters
input = "/home/daniel/mpi/conference/ohbm/ohbm_2019/poster/odc_exp/data/rh.spmT_left_right_GE_EPI5_def_layer9_sigma0_grid.nii"
x_grid = 0.25 # grid resolution in mm (x-direction)
y_grid = 0.25 # grid resolution in mm (y-direction)
name_output = "GE_EPI5" # basename of output
path_output = "/home/daniel/mpi/conference/ohbm/ohbm_2019/poster/odc_exp/results" # path where output is saved

""" do not edit below """

# make output folder
if not os.path.exists(path_output):
    os.mkdir(path_output)

# get number of layers
data = nb.load(input).get_fdata()
nlayer = np.shape(data)[2]

# hemisphere
hemi = os.path.splitext(os.path.splitext(os.path.basename(input))[0])[0]

# considered rotation angles for generation of projection lines
phi = 10 * np.arange(36)

k_fft_phi = np.zeros((nlayer, len(phi)))
P_fft_phi = np.zeros((nlayer, len(phi)))
d_acorr_phi = np.zeros((nlayer, len(phi)))
P_acorr_phi = np.zeros((nlayer, len(phi)))
fwhm_acorr_phi = np.zeros((nlayer, len(phi)))
x_fft_0 = []
y_fft_0 = []
x_fft_90 = []
y_fft_90 = []
x_acorr_0 = []
y_acorr_0 = []
x_acorr_90 = []
y_acorr_90 = []
for i in range(nlayer):
    
    # load nifti
    data = nb.load(input).get_fdata()[:,:,i]

    # field of view of input data in mm
    FOVx = np.shape(data)[0] * x_grid
    FOVy = np.shape(data)[1] * y_grid

    # get pca
    _, _, x_minor, y_minor = get_pca(data)
    
    # compute for each rotation angle
    for j in range(len(phi)):
    
        # rotate grid
        x_temp = x_minor*np.cos(phi[j] / 180 * np.pi) - y_minor*np.sin(phi[j] / 180 * np.pi)
        y_temp = x_minor*np.sin(phi[j] / 180 * np.pi) + y_minor*np.cos(phi[j] / 180 * np.pi)
    
        # analyze fourier spectrum
        k_fft, P_fft, x_fft, y_fft = analyze_fft(data, FOVx, FOVy, x_temp, y_temp)
    
        # analyze autocorrelation
        fwhm_acorr, d_acorr, P_acorr, x_acorr, y_acorr = analyze_acorr(data, FOVx, FOVy, x_temp, y_temp)
        
        # list result
        k_fft_phi[i,j] = float(k_fft)
        P_fft_phi[i,j] = float(P_fft)
        d_acorr_phi[i,j] = float(d_acorr)
        P_acorr_phi[i,j] = float(P_acorr)
        fwhm_acorr_phi[i,j] = float(fwhm_acorr)
        
        # get example data for minor and major axes
        if phi[j] == 0:
            x_fft_0.append(x_fft)
            y_fft_0.append(y_fft)
            x_acorr_0.append(x_acorr)
            y_acorr_0.append(y_acorr)
        
        if phi[j] == 90:
            x_fft_90.append(x_fft)
            y_fft_90.append(y_fft)
            x_acorr_90.append(x_acorr)
            y_acorr_90.append(y_acorr)

# save variables
np.savez(os.path.join(path_output,hemi+"."+name_output),
         x_fft_0=x_fft_0, x_fft_90=x_fft_90, x_acorr_0=x_acorr_0, x_acorr_90=x_acorr_90,
         y_fft_0=y_fft_0, y_fft_90=y_fft_90, y_acorr_0=y_acorr_0, y_acorr_90=y_acorr_90,
         k_fft_phi=k_fft_phi, P_fft_phi=P_fft_phi, 
         d_acorr_phi=d_acorr_phi, P_acorr_phi=P_acorr_phi, fwhm_acorr_phi=fwhm_acorr_phi)