# -*- coding: utf-8 -*-
"""
Simulation of imaging MTF

The purpose of the following script is to simulate the imaging modulation
transfer function from a GE-EPI from imaging of a thin slid. The slid is
parallel to the x-direction (read direction). The following steps are done:
    1. generate input image
    2. simulate MR sampling including k-space truncation, pF (zero filling) and GRAPPA and T2* decay
    3. zero pad raw data matrix
    3. get magnitude image
    4. plot projection of mtf (included is correction for slid width)
    5. plot input and output images

"""

# python standard library inputs
import os

# external inputs
import numpy as np
import nibabel as nb
from numpy.fft import fftshift, fft2, ifft2
import matplotlib.pyplot as plt
from matplotlib import rc

# parameters
Nx_sim = 1024
Ny_sim = 1024
Nx_mri = 184
Ny_mri = 184
Nx_pad = 512
Ny_pad = 512
FOVx = 148  # in mm
FOVy = 148  # in mm
slot_x = 70  # in mm
slot_y = 0.5  # in mm

T2star = 27  # in ms
esp = 1  # in ms
TE = 24  # in ms
grappa = 3
pf = 6 / 8

path_output = "/home/daniel/Schreibtisch"

# do not edit below

# font parameters for plots
rc('font', **{'family': 'serif', 'serif': ['Palatino']})
rc('text', usetex=True)

# generate image
arr = np.zeros((Ny_sim, Nx_sim))

slot_length = np.round(slot_x * Nx_sim / FOVx)
slot_width = np.round(slot_y * Ny_sim / FOVy)

x1 = np.round(Nx_sim / 2).astype(int) - np.round(slot_length / 2).astype(int)
x2 = np.round(x1 + slot_length).astype(int)
y1 = np.round(Ny_sim / 2).astype(int) - np.round(slot_width / 2).astype(int)
y2 = np.round(y1 + slot_width).astype(int)

arr[y1:y2, x1:x2] = 1

# FFT
y_fft = fft2(arr)

# truncate k-space
kx_sample = np.round(Nx_mri / 2).astype(int)
ky_sample = np.round(Ny_mri / 2).astype(int)

ymri_fft = np.zeros((Ny_mri, Nx_mri), dtype=complex)
ymri_fft[:ky_sample, :kx_sample] = y_fft[:ky_sample, :kx_sample]
ymri_fft[:ky_sample, -1:-kx_sample - 1:-1] = y_fft[:ky_sample, -1:-kx_sample - 1:-1]
ymri_fft[-1:-ky_sample - 1:-1, :kx_sample] = y_fft[-1:-ky_sample - 1:-1, :kx_sample]
ymri_fft[-1:-ky_sample - 1:-1, -1:-kx_sample - 1:-1] = y_fft[-1:-ky_sample - 1:-1, -1:-kx_sample - 1:-1]
ymri_fft = fftshift(ymri_fft)  # higher frequencies to edge

# k-space weighting
Ny_grappa = np.ceil(Ny_mri / grappa).astype(int)
t = np.arange(0, Ny_grappa, esp) - np.floor(Ny_grappa / 2).astype(int) + TE
ymri_accel = np.zeros((Ny_grappa, Nx_mri))
for j in range(Ny_grappa):
    ymri_accel[j, :] = np.exp(-t[j] / T2star)
ymri_accel = np.repeat(ymri_accel, grappa, axis=0)
ymri_accel = ymri_accel[:np.shape(ymri_fft)[0], :np.shape(ymri_fft)[1]]
ymri_accel[:np.round((1 - pf) * Ny_mri).astype(int), :] = 0  # partial Fourier
ymri_fft = ymri_fft * ymri_accel

# pad data matrix
pad_x = int(Nx_pad / 2 - Nx_mri / 2)
pad_y = int(Ny_pad / 2 - Ny_mri / 2)
ymri_fft = np.pad(ymri_fft, (pad_y, pad_x), mode="constant")

# iFFT
ymri_fft = fftshift(ymri_fft)  # lower frequencies to edge
ymri_magn = np.abs(ifft2(ymri_fft))

# plot normalized mtf
fig, ax = plt.subplots()
f = fftshift(np.linspace(-Ny_pad / FOVy, Ny_pad / FOVy, Ny_pad))
mtf = np.abs(fft2(ymri_magn))
f_ny = np.sqrt(Nx_mri ** 2 + Ny_mri ** 2) / (2 * np.sqrt(FOVx ** 2 + FOVy ** 2))  # nyquist frequency
x = f[:np.round(Ny_pad / 2).astype(int)]
y = mtf[:np.round(Ny_pad / 2).astype(int), 0]
y = y / np.nanmax(y)
y = y[x < 1] * np.pi * x[x < 1] * slot_y / np.sin(np.pi * x[x < 1] * slot_y)
x = x[x < 1]
ax.plot(x, y, linewidth=1.0)
ax.set_xlabel("Spatial frequency in cycles/mm")
ax.set_ylabel("Normalized MTF in a.u.")
ax.axvline(x=f_ny, ymin=0, ymax=1, linestyle='--', color='red')
ax.text(f_ny - 0.05, 1.08, r"$f_{\mathrm{Nyquist}}$")
ax.legend(["pF = 6/8, iPAT = 3"], frameon=False)
fig.savefig(os.path.join(path_output, "mtf_example.svg"), format='svg', bbox_inches='tight')
plt.show()

# save input and output images as nifti
output = nb.Nifti1Image(arr, np.eye(4), nb.Nifti1Header())
nb.save(output, os.path.join(path_output, "input.nii"))

output = nb.Nifti1Image(ymri_magn, np.eye(4), nb.Nifti1Header())
nb.save(output, os.path.join(path_output, "output.nii"))
