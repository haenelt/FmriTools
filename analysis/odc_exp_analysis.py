"""
Analysis of simulated ocular dominance columns

The purpose of the following script is to simultae ocular dominance columns and estimate their 
column and column spacing. To do so, the power spectrum and the normalized autocorrelation (NAC) map 
are computed and evaluated along different projection lines. Especially, principal axes are first 
estimated from the thresholded Fourier spectogram and width and spacing results are saved for these
directions. In total, the following quantities are computed:

(1) column width is estimated by computing the FWHM of the NAC central peak.
(2) column spacing is estimated by computing the distance from the NAC central peak to its first 
neighbor.
(3) take amplitude of first neighbor peak as repetition marker?
(3) column spacing is further estimated by computing the spatial frequency of the maximum Fourier
component.
(4) take the maximum power (relative to central peak) as repetition marker?

Line and polar plots are saved from all measures. Further below, not all parameters are needed for
this simulation since we only use the neural map, i.e., no MRI parameters (Nx_mri, Ny_mri),
bold parameters (beta, fwhm_bold, fwhm_noise) or occlusion parameters (a, b, theta) are needed.

created by Daniel Haenelt
Date created: 15-04-2019
Last modified: 15-04-2019
"""
import os
import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
from lib.analysis import get_pca
from lib.analysis import analyze_fft
from lib.analysis import analyze_acorr

# parameters for ODC analysis
input = "/home/daniel/Schreibtisch/intermediate/img/odc_exp_analysis/data/lh.spmT_left_right_GE_EPI3_def_layer9_sigma0_grid.nii"
phi = [0,90] # considered angles for generation of projection lines
path_output = "/home/daniel/Schreibtisch/test" # path where output is saved

FOVx_mri = 148
FOVy_mri = 148
Nx_mri = 184
Ny_mri = 184

""" do not edit below """

# make output folder
if not os.path.exists(path_output):
    os.mkdir(path_output)

# nyquist frequency
f_ny = np.sqrt(Nx_mri**2+Ny_mri**2) / np.sqrt(FOVx_mri**2+FOVy_mri**2)

# load nifti
data = nb.load(input).get_fdata()[:,:,0]

FOVx = np.shape(data)[0]*0.25
FOVy = np.shape(data)[1]*0.25

k_fft_res = np.zeros(len(phi))
P_fft_res = np.zeros(len(phi))
d_acorr_res= np.zeros(len(phi))
P_acorr_res = np.zeros(len(phi))
fwhm_acorr_res = np.zeros(len(phi))

# get pca
_, _, x_minor, y_minor = get_pca(data)
    
for j in range(len(phi)):
    
    # rotate grid
    x_temp = x_minor*np.cos(phi[j] / 180 * np.pi) - y_minor*np.sin(phi[j] / 180 * np.pi)
    y_temp = x_minor*np.sin(phi[j] / 180 * np.pi) + y_minor*np.cos(phi[j] / 180 * np.pi)
    
    # analyze fourier spectrum
    k_fft, P_fft, _, _ = analyze_fft(data, FOVx, FOVy, x_temp, y_temp, f_cut=0.05)
    
    # analyze autocorrelation
    fwhm_acorr, d_acorr, P_acorr, _, _ = analyze_acorr(data, FOVx, FOVy, x_temp, y_temp, 
                                                       0, None)
        
    # list result
    k_fft_res[j] = k_fft
    P_fft_res[j] = P_fft
    d_acorr_res[j] = d_acorr
    P_acorr_res[j] = P_acorr
    fwhm_acorr_res[j] = fwhm_acorr

"""
Example plots
"""

# get pca
x_major, y_major, x_minor, y_minor = get_pca(data)

# analyze fft
_, _, x1, y1 = analyze_fft(data, FOVx, FOVy, x_major, y_major, 10, None)
_, _, x2, y2 = analyze_fft(data, FOVx, FOVy, x_minor, y_minor, 10, None)

# fft
y1 = y1[x1 < 2]
x1 = x1[x1 < 2]
y2 = y2[x2 < 2]
x2 = x2[x2 < 2]

fig, ax = plt.subplots()
ax.plot(x1, y1)
ax.plot(x2, y2)
ax.axvline(x=1/f_ny, ymin=0, ymax=1)
ax.set_xlabel("Spatial frequency in cycles/mm")
ax.set_ylabel("Power spectrum in a.u.")
ax.legend(["major axis","minor axis"])
fig.savefig(os.path.join(path_output,"fft_example.png"), bbox_inches="tight")
plt.show()

# analyze autocorrelation
_, _, _, x1, y1 = analyze_acorr(data, FOVx, FOVy, x_major, y_major, 0.01, None)
_, _, _, x2, y2 = analyze_acorr(data, FOVx, FOVy, x_minor, y_minor, 0.01, None)

# nac
fig, ax = plt.subplots()
ax.plot(x1, y1)
ax.plot(x2, y2)
ax.set_xlabel("Lag in mm")
ax.set_ylabel("NAC")
ax.legend(["major axis","minor axis"])
fig.savefig(os.path.join(path_output,"acorr_example.png"), bbox_inches="tight")
plt.show()

print(k_fft_res)
