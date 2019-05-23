"""
Analysis of simulated ocular dominance columns

The purpose of the following script is to simulate ocular dominance columns and estimate their 
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

Further below, not all parameters are needed for this simulation since we only use the neural map, 
i.e., no MRI parameters (Nx_mri, Ny_mri), bold parameters (beta, fwhm_bold, fwhm_noise) or occlusion 
parameters (a_mask, b_mask, theta_mask) are needed.

created by Daniel Haenelt
Date created: 14-04-2019
Last modified: 22-05-2019
"""
import os
import numpy as np
from lib.simulation.odc import odc_2d
from lib.analysis import get_pca
from lib.analysis import analyze_fft
from lib.analysis import analyze_acorr

# parameters for ODC image matrix
Nx_sim = 1024
Ny_sim = 1024
FOVx = 20
FOVy = 20
Nx_mri = 25
Ny_mri = 25

# parameters for ODC shape
rho = 0.5
delta = 0.3
epsilon = 0.4
theta = 90
alpha = 4

# parameters for bold response
beta = 1
fwhm_bold = 1.02
fwhm_noise = 0

# parameters for edge occlusion
a_mask = 1000
b_mask = 1000
theta_mask = 0

# parameters for ODC analysis
niter = 1000 # number of iterations
name_output = "sim" # basename of output
path_output = "/data/pt_01880/odc_sim" # path where output is saved

""" do not edit below """

# make output folder
if not os.path.exists(path_output):
    os.mkdir(path_output)

# considered rotation angles for generation of projection lines
phi = 10 * np.arange(36)

k_fft_phi = np.zeros((niter, len(phi)))
P_fft_phi = np.zeros((niter, len(phi)))
d_acorr_phi= np.zeros((niter, len(phi)))
P_acorr_phi = np.zeros((niter, len(phi)))
fwhm_acorr_phi = np.zeros((niter, len(phi)))
x_fft_0 = []
y_fft_0 = []
x_fft_90 = []
y_fft_90 = []
x_acorr_0 = []
y_acorr_0 = []
x_acorr_90 = []
y_acorr_90 = []
for i in range(niter):
    
    # get odc pattern
    _, neural, _, _, _, _ = odc_2d(Nx_sim, Ny_sim, FOVx, FOVy, Nx_mri, Ny_mri, rho, delta, epsilon, 
                                   theta, alpha, beta, fwhm_bold, fwhm_noise, a_mask, b_mask, 
                                   theta_mask, False)

    # get pca
    _, _, x_minor, y_minor = get_pca(neural)

    for j in range(len(phi)):

        # rotate grid
        x_temp = x_minor*np.cos(phi[j] / 180 * np.pi) - y_minor*np.sin(phi[j] / 180 * np.pi)
        y_temp = x_minor*np.sin(phi[j] / 180 * np.pi) + y_minor*np.cos(phi[j] / 180 * np.pi)
   
        # analyze fourier spectrum
        k_fft, P_fft, x_fft, y_fft = analyze_fft(neural, FOVx, FOVy, x_temp, y_temp)
    
        # analyze autocorrelation
        fwhm_acorr, d_acorr, P_acorr, x_acorr, y_acorr = analyze_acorr(neural, FOVx, FOVy, x_temp, y_temp)
        
        # list result
        k_fft_phi[i,j] = k_fft
        P_fft_phi[i,j] = P_fft
        d_acorr_phi[i,j] = d_acorr
        P_acorr_phi[i,j] = P_acorr
        fwhm_acorr_phi[i,j] = fwhm_acorr
        
        # get example data for minor and major axes
        if i == 0 and phi[j] == 0:
            x_fft_0.append(x_fft)
            y_fft_0.append(y_fft)
            x_acorr_0.append(x_acorr)
            y_acorr_0.append(y_acorr)
        
        if i == 0 and phi[j] == 90:
            x_fft_90.append(x_fft)
            y_fft_90.append(y_fft)
            x_acorr_90.append(x_acorr)
            y_acorr_90.append(y_acorr)
    
# mean across iterations
k_fft_phi_mean = np.nanmean(k_fft_phi,0)
P_fft_phi_mean = np.nanmean(P_fft_phi,0)
d_acorr_phi_mean = np.nanmean(d_acorr_phi,0)
P_acorr_phi_mean = np.nanmean(P_acorr_phi,0)
fwhm_acorr_phi_mean = np.nanmean(fwhm_acorr_phi,0)
    
# std across iterations
k_fft_phi_std = np.nanstd(k_fft_phi,0)
P_fft_phi_std = np.nanstd(P_fft_phi,0)
d_acorr_phi_std = np.nanstd(d_acorr_phi,0)
P_acorr_phi_std = np.nanstd(P_acorr_phi,0)
fwhm_acorr_phi_std = np.nanstd(fwhm_acorr_phi,0)

# save variables
np.savez(os.path.join(path_output,name_output),
         x_fft_0=x_fft_0, x_fft_90=x_fft_90, x_acorr_0=x_acorr_0, x_acorr_90=x_acorr_90,
         y_fft_0=y_fft_0, y_fft_90=y_fft_90, y_acorr_0=y_acorr_0, y_acorr_90=y_acorr_90,
         k_fft_phi_mean=k_fft_phi_mean, P_fft_phi_mean=P_fft_phi_mean,
         d_acorr_phi_mean=d_acorr_phi_mean, P_acorr_phi_mean=P_acorr_phi_mean, fwhm_acorr_phi_mean=fwhm_acorr_phi_mean,
         k_fft_phi_std=k_fft_phi_std, P_fft_phi_std=P_fft_phi_std,
         d_acorr_phi_std=d_acorr_phi_std, P_acorr_phi_std=P_acorr_phi_std, fwhm_acorr_phi_std=fwhm_acorr_phi_std)