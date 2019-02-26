"""
Generate 2D ODC

This script creates a two-dimensional ODC pattern.

created by Daniel Haenelt
Date created: 18-02-2019
Last modified: 18-02-2019
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate
from numpy.fft import fft2, fftshift
from lib.simulation.odc import odc_2d

"""
2D parameters
"""

Nx_sim = 1024
Ny_sim = 1024
FOVx = 20
FOVy = 20
Nx_mri = 25
Ny_mri = 25

rho = 0.5
delta = 0.3
epsilon = 0.4
theta = 90
alpha = 4

beta = 1
fwhm_bold = 1.02
fwhm_noise = 0

"""
Generate pattern
"""

_, neural, bold, mri, _, _ = odc_2d(Nx_sim, Ny_sim, FOVx, FOVy, Nx_mri, Ny_mri, rho, delta, epsilon, 
                                    theta, alpha, beta, fwhm_bold, fwhm_noise, False)

"""
Plot
"""

# plot neural map
_, ax = plt.subplots()
ax.imshow(neural, extent=[0,FOVx,0,FOVy])
ax.set_xlabel("Distance in mm")
ax.set_ylabel("Distance in mm")
ax.set_title("Simulated pattern")

# plot sampled pattern
_, ax = plt.subplots()
ax.imshow(mri,extent=[0,FOVx,0,FOVy])
ax.set_xlabel("Distance in mm")
ax.set_ylabel("Distance in mm")
ax.set_title("Sampled pattern")

# autocorrelation
mri_auto = correlate(mri, mri, "same")
_, ax = plt.subplots()
ax.imshow(mri_auto, extent=[-FOVx/2,FOVx/2,-FOVy/2,FOVy/2])
ax.set_xlabel("Lag in mm")
ax.set_ylabel("Lag in mm")
ax.set_title("Autocorrelation")

# fft
mri_fft = fftshift(np.abs(fft2(mri))**2)
_, ax = plt.subplots()
ticks = [-Nx_mri/(2*FOVx),Nx_mri/(2*FOVx),-Ny_mri/(2*FOVy),Ny_mri/(2*FOVy)]
if not Nx_mri % 2:
    ticks[0] = (-Nx_mri-1)/(2*FOVx)
    ticks[1] = (Nx_mri-1)/(2*FOVy)
if not Ny_mri % 2:
    ticks[2] = (-Ny_mri-1)/(2*FOVy)
    ticks[3] = (Ny_mri-1)/(2*FOVy)
ax.imshow(mri_fft, extent=ticks)
ax.set_xlabel("Spatial frequency in cycles/mm")
ax.set_ylabel("Spatial frequency in cycles/mm")
ax.set_title("Spatial frequency representation")