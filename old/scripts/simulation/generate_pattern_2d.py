# -*- coding: utf-8 -*-
"""
Generate 2D pattern

This script creates a two-dimensional sinusoidal or rectangular pattern.

"""

import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import fft2, fftshift
from scipy.signal import correlate

from ..utils.simulation import pattern_2d

# 2D parameters
Nx_sim = 1024
Ny_sim = 1024
FOVx = 20
FOVy = 20
Nx_mri = 25
Ny_mri = 25

omega_x = 0.1
omega_y = 0
phi_x = 0
phi_y = 0
theta = 0
rect_shape = True

beta = 1
fwhm_bold = 0
fwhm_noise = 0

a_mask = 1000
b_mask = 1000
theta_mask = 0

# do not edit below

# Generate pattern
neural, bold, mri, _ = pattern_2d(
    Nx_sim,
    Ny_sim,
    FOVx,
    FOVy,
    Nx_mri,
    Ny_mri,
    omega_x,
    omega_y,
    phi_x,
    phi_y,
    theta,
    rect_shape,
    beta,
    fwhm_bold,
    fwhm_noise,
    a_mask,
    b_mask,
    theta_mask,
)

# Plot
# plot neural map
_, ax = plt.subplots()
ax.imshow(neural, extent=[0, FOVx, 0, FOVy])
ax.set_xlabel("Distance in mm")
ax.set_ylabel("Distance in mm")
ax.set_title("Simulated pattern")

# plot sampled pattern
_, ax = plt.subplots()
ax.imshow(mri, extent=[0, FOVx, 0, FOVy])
ax.set_xlabel("Distance in mm")
ax.set_ylabel("Distance in mm")
ax.set_title("Sampled pattern")

# autocorrelation
mri_auto = correlate(mri, mri, "same")
_, ax = plt.subplots()
ax.imshow(mri_auto, extent=[-FOVx / 2, FOVx / 2, -FOVy / 2, FOVy / 2])
ax.set_xlabel("Lag in mm")
ax.set_ylabel("Lag in mm")
ax.set_title("Autocorrelation")

# fft
mri_fft = fftshift(np.abs(fft2(mri)) ** 2)
_, ax = plt.subplots()
ticks = [
    -Nx_mri / (2 * FOVx),
    Nx_mri / (2 * FOVx),
    -Ny_mri / (2 * FOVy),
    Ny_mri / (2 * FOVy),
]
if not Nx_mri % 2:
    ticks[0] = (-Nx_mri - 1) / (2 * FOVx)
    ticks[1] = (Nx_mri - 1) / (2 * FOVy)
if not Ny_mri % 2:
    ticks[2] = (-Ny_mri - 1) / (2 * FOVy)
    ticks[3] = (Ny_mri - 1) / (2 * FOVy)
ax.imshow(mri_fft, extent=ticks)
ax.set_xlabel("Spatial frequency in cycles/mm")
ax.set_ylabel("Spatial frequency in cycles/mm")
ax.set_title("Spatial frequency representation")
