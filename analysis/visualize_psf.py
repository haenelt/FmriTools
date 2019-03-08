"""
Visualization of PSF data across cortical depth.

created by Daniel Haenelt
Date created: 08-03-2019
Last modified: 08-03-2019
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

# input
path = "/data/pt_01880/test" # path to saved input values
nlayer = 10 # number of layers
apply_gaussian = True

""" do not edit below """

# initialize array
x = np.load(os.path.join(path,"x_0.npy"))
y = np.load(os.path.join(path,"y_0.npy"))

# get array
X = x
Y = np.arange(nlayer-1,-1,-1)
X, Y = np.meshgrid(X,Y)

# get data
Z = np.zeros((nlayer,len(y)))
for i in range(nlayer):
    Z[i,:] = 100 * np.load(os.path.join(path,"y_"+str(i)+".npy"))

if apply_gaussian:
    Z = gaussian_filter(Z, 1.0) # with sigma 1.0

# figure 1 (2D plot of MTf data)
fig = plt.figure()
im = plt.pcolor(X, Y, Z,cmap='hot')
cbar = plt.colorbar(im) # adding the colobar on the right
cbar.set_label("Coherence value in percent")
plt.xlabel("Cortical frequency in cycles/mm")
plt.ylabel("Cortical depth (WM -> CSF)")
plt.contour(X, Y, Z, 10, colors='g',interpolation='none')
plt.savefig(os.path.join(path,"normal_mtf.png"),dpi=100)
plt.show()

# figure 2 (1D plot of estimated PSF)
sigma = np.zeros(nlayer)
for i in range(nlayer):
    sigma[i] = np.load(os.path.join(path,"sigma_"+str(i)+".npy"))

fwhm = 1/np.abs(2*np.sqrt(2*np.log(2))*sigma)

fig2 = plt.figure()
plt.plot(fwhm[-1:0:-1])
plt.xlabel("Cortical depth (WM -> CSF)")
plt.ylabel("FWHM in mm")
plt.savefig(os.path.join(path,"normal_psf.png"))