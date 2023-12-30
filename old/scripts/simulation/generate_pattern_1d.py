# -*- coding: utf-8 -*-
"""
Generate 1D pattern

This script creates a one-dimensional sinusoidal or rectangular pattern.

"""

import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import fft
from scipy.signal import correlate

from ..simulation.pattern import pattern_1d

# 1D parameters
N_sim = 1024
FOV = 20
N_mri = 500

omega = 0.5
phi = 0
rect_shape = True

beta = 1
fwhm_bold = 0
fwhm_noise = 0

a = 1000
b = 1000

# do not edit below

# Generate pattern
neural, bold, mri, _ = pattern_1d(
    N_sim, FOV, N_mri, omega, phi, rect_shape, beta, fwhm_bold, fwhm_noise, a, b
)

# Plot
# plot neural map
x = np.linspace(0, 1, N_sim) * FOV
_, ax = plt.subplots()
ax.plot(x, neural)
ax.set_xlabel("Distance in mm")
ax.set_ylabel("Signal change in a.u.")
ax.set_title("Simulated pattern")

# plot sampled pattern
x = np.linspace(0, 1, N_mri) * FOV
_, ax = plt.subplots()
ax.plot(x, mri)
ax.set_xlabel("Distance in mm")
ax.set_ylabel("Signal change in a.u.")
ax.set_title("Sampled pattern")

# autocorrelation
x = np.linspace(-1, 1, N_mri) * FOV
mri_auto = correlate(mri, mri, "same")
_, ax = plt.subplots()
ax.plot(x, mri_auto)
ax.set_xlabel("Lag in mm")
ax.set_ylabel("ACF")
ax.set_title("Autocorrelation")

# fft
x = np.arange(0, np.floor(N_mri / 2)) / FOV
mri_fft = np.abs(fft(mri)) ** 2
_, ax = plt.subplots()
ax.plot(x, mri_fft[: np.floor(N_mri / 2).astype(int)])
ax.set_xlabel("Spatial frequency in cycles/mm")
ax.set_ylabel("FFT")
ax.set_title("Spatial frequency representation")
