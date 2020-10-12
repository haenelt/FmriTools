# -*- coding: utf-8 -*-

# python standard library inputs
import os

# external inputs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


"""
Plots for autocorrelation analysis (simulation)

The purpose of the following script is to generate exemplary plots for the 
autocorrelation analysis of ODC data (averaged over 1000 iterations). The 
following plots are produced:
    1. phi vs. fwhm
    2. phi vs. k_fft 
    3. phi vs. P_fft
    4. phi vs. d_acorr
    5. phi vs. P_acorr
    6. fft (for single iteration)
    7. autocorrelation (for single iteration)

created by Daniel Haenelt
Date created: 18-04-2019
Last modified: 12-10-2020
"""

# load data from simulation
Y = np.load("/home/daniel/mpi/conference/ohbm/ohbm_2019/poster/odc_sim/sim.npz")

path_output = "/home/daniel/mpi/conference/ohbm/ohbm_2019/poster/odc_sim"

# do not edit below

# font parameters for plots
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

# phi vs. fwhm, k_fft, P_fft, d_acorr, P_acorr

# fwhm
phi = 10*np.arange(36)
y_mean = Y["fwhm_acorr_phi_mean"]
y_std = Y["fwhm_acorr_phi_std"]


fig, ax = plt.subplots()
ax.plot(phi, y_mean, "r")
ax.fill_between(phi, y_mean-y_std, y_mean+y_std, alpha=0.2, facecolor='b')
ax.set_xlabel("Rotation angle in deg")
ax.set_ylabel("Column width in mm")
ax.set_title("ODC along different projections")
fig.savefig(os.path.join(path_output,"sim_fwhm_vs_phi.svg"), format='svg', bbox_inches='tight')
plt.show()

# k_fft
phi = 10*np.arange(36)
y_mean = Y["k_fft_phi_mean"]
y_std = Y["k_fft_phi_std"]

fig, ax = plt.subplots()
ax.plot(phi, y_mean, "r")
ax.fill_between(phi, y_mean-y_std, y_mean+y_std, alpha=0.2, facecolor='b')
ax.set_xlabel("Rotation angle in deg")
ax.set_ylabel("Peak spatial frequency in cycles/mm")
ax.set_title("ODCs along different projections")
fig.savefig(os.path.join(path_output,"sim_kfft_vs_phi.svg"), format='svg', bbox_inches='tight')
plt.show()

# P_fft
phi = 10*np.arange(36)
y_mean = Y["P_fft_phi_mean"]
y_std = Y["P_fft_phi_std"]

fig, ax = plt.subplots()
ax.plot(phi, y_mean, "r")
ax.fill_between(phi, y_mean-y_std, y_mean+y_std, alpha=0.2, facecolor='b')
ax.set_ylim([0,20000])
ax.set_xlabel("Rotation angle in deg")
ax.set_ylabel("Peak spatial frequency magnitude in a.u.")
ax.set_title("ODCs along different projections")
#fig.savefig(os.path.join(path_output,"sim_Pfft_vs_phi.svg"), format='svg', bbox_inches='tight')
fig.savefig(os.path.join(path_output,"sim_Pfft_vs_phi_threshold.svg"), format='svg', bbox_inches='tight')
plt.show()

# d_acorr
phi = 10*np.arange(36)
y_mean = Y["d_acorr_phi_mean"]
y_std = Y["d_acorr_phi_std"]

fig, ax = plt.subplots()
ax.plot(phi, y_mean, "r")
ax.fill_between(phi, y_mean-y_std, y_mean+y_std, alpha=0.2, facecolor='b')
ax.set_xlabel("Rotation angle in deg")
ax.set_ylabel("Nearest neighbor in mm")
ax.set_title("ODCs along different projections")
fig.savefig(os.path.join(path_output,"sim_dacorr_vs_phi.svg"), format='svg', bbox_inches='tight')
plt.show()

# P_acorr
phi = 10*np.arange(36)
y_mean = Y["P_acorr_phi_mean"]
y_std = Y["P_acorr_phi_std"]

fig, ax = plt.subplots()
ax.plot(phi, y_mean, "r")
ax.fill_between(phi, y_mean-y_std, y_mean+y_std, alpha=0.2, facecolor='b')
ax.set_xlabel("Rotation angle in deg")
ax.set_ylabel("Nearest neighbor magnitude in a.u.")
ax.set_title("ODCs along different projections")
fig.savefig(os.path.join(path_output,"sim_Pacorr_vs_phi.svg"), format='svg', bbox_inches='tight')
plt.show()

# FFT
y1 = Y["y_fft_0"]
x1 = Y["x_fft_0"]
y1 = y1[x1 < 2]
x1 = x1[x1 < 2]

y2 = Y["y_fft_90"]
x2 = Y["x_fft_90"]
y2 = y2[x2 < 2]
x2 = x2[x2 < 2]

# fwhm plot
fig, ax = plt.subplots()
ax.plot(x1,y1,linewidth=1.0)
ax.plot(x2,y2,linewidth=1.0)
ax.set_xlabel("Spatial frequency in cycles/mm")
ax.set_ylabel("Normalized magnitude in a.u.")
ax.set_title("Fourier spectrum")
ax.legend(["minor axis", "major axis"], frameon=False)
fig.savefig(os.path.join(path_output,"sim_fft_example.svg"), format='svg', bbox_inches='tight')
plt.show()

# autocorrelation
x1 = Y["x_acorr_0"][0]
y1 = Y["y_acorr_0"][0]
x2 = Y["x_acorr_90"][0]
y2 = Y["y_acorr_90"][0]

# fwhm plot
fig, ax = plt.subplots()
ax.plot(x1,y1,linewidth=1.0)
ax.plot(x2,y2,linewidth=1.0)
ax.set_xlabel("Lag in mm")
ax.set_ylabel("NAC in a.u.")
ax.set_title("Autocorrelation")
ax.legend(["minor axis", "major axis"], frameon=False)
fig.savefig(os.path.join(path_output,"sim_acorr_example.svg"), format='svg', bbox_inches='tight')
plt.show()
