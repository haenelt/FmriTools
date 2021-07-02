# -*- coding: utf-8 -*-
"""
Visualization of PSF data across cortical depth.

"""

# python standard library inputs
import os

# external inputs
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

# input
path = "/data/pt_01880/test"  # path to saved input values
nlayer = 10  # number of layers
apply_gaussian = True
gaussian = True
nsigma = 1

# do not edit below

# get array
X = np.load(os.path.join(path, "F_out.npy"))
Y = np.arange(nlayer - 1, -1, -1)
X, Y = np.meshgrid(X, Y)

# get data
Z = 100 * np.load(os.path.join(path, "M_out.npy"))

# apply gaussian filter
if apply_gaussian:
    Z = gaussian_filter(Z, 1.0)  # with sigma 1.0

# figure 1 (2D plot of MTf data)
fig = plt.figure()
plt.rcParams["font.size"] = 10
im = plt.pcolor(X, Y, Z, cmap='hot')
cbar = plt.colorbar(im)  # adding the colobar on the right
cbar.set_label("Coherence value in percent")
plt.xlabel("Cortical frequency in cycles/mm")
plt.ylabel("Cortical depth (WM -> CSF)")
plt.contour(X, Y, Z, 10, colors='g', interpolation='none')
plt.savefig(os.path.join(path, "normal_mtf.png"), dpi=100)
plt.show()

# figure 2 (1D plot of estimated PSF)
sigma = np.load(os.path.join(path, "sigma_out.npy"))
sigma = sigma[::-1]  # reverse to get have array in WM -> CSF direction

if gaussian:
    fwhm = 1 / np.abs(2 * np.sqrt(2 * np.log(2)) * sigma[:, 0])
    fwhm_dw = 1 / np.abs(2 * np.sqrt(2 * np.log(2)) * sigma[:, 1])
    fwhm_up = 1 / np.abs(2 * np.sqrt(2 * np.log(2)) * sigma[:, 2])
else:
    fwhm = 2 * sigma[:, 0]
    fwhm_dw = 2 * sigma[:, 1]
    fwhm_up = 2 * sigma[:, 2]

fig, ax = plt.subplots(1)
plt.rcParams["font.size"] = 10
plt.plot(fwhm, "r", lw=1)
ax.fill_between(np.arange(0, 10), fwhm_up, fwhm_dw, alpha=.25,
                label=str(nsigma) + "-sigma interval")
plt.xlabel("Cortical depth (WM -> CSF)")
plt.ylabel("FWHM in mm")
plt.savefig(os.path.join(path, "normal_psf.png"))
plt.legend(loc="upper left", fontsize=12, frameon=False)
