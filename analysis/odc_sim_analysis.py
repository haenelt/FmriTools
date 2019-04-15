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
Date created: 14-04-2019
Last modified: 14-04-2019
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import rotate
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
a = 1000
b = 1000
theta = 0

# parameters for ODC analysis
phi = 10*np.arange(36) # considered angles for generation of projection lines
niter = 2 # number of iterations
path_output = "/home/daniel/Schreibtisch/test" # path where output is saved

""" do not edit below """

# make output folder
if not os.path.exists(path_output):
    os.mkdir(path_output)

k_fft_res = np.zeros((niter,len(phi)))
P_fft_res = np.zeros((niter,len(phi)))
d_acorr_res= np.zeros((niter,len(phi)))
P_acorr_res = np.zeros((niter,len(phi)))
fwhm_acorr_res = np.zeros((niter,len(phi)))
for i in range(niter):
    
    # get odc pattern
    _, neural, _, _, _, _ = odc_2d(Nx_sim, Ny_sim, FOVx, FOVy, Nx_mri, Ny_mri, rho, delta, epsilon, 
                                   theta, alpha, beta, fwhm_bold, fwhm_noise, a, b, theta, False)

    # get pca
    _, _, x_minor, y_minor = get_pca(neural)

    # rotate grid
    for j in range(len(phi)):
        neural_temp = rotate(neural,phi[j])
        if np.mod(phi[j],90) != 0:
            x_cut = int(np.shape(neural_temp)[0]/2-np.shape(neural)[0]/2)
            y_cut = int(np.shape(neural_temp)[1]/2-np.shape(neural)[1]/2)
            neural_temp = neural_temp[x_cut:-x_cut,y_cut:-y_cut]
    
        # analyze fourier spectrum
        k_fft, P_fft, _, _ = analyze_fft(neural_temp, FOVx, FOVy, x_minor, y_minor, 10, None)
    
        # analyze autocorrelation
        fwhm_acorr, d_acorr, P_acorr, _, _ = analyze_acorr(neural_temp, FOVx, FOVy, x_minor, 
                                                           y_minor, 0.01, None)
        
        # list result
        k_fft_res[i,j] = k_fft
        P_fft_res[i,j] = P_fft
        d_acorr_res[i,j] = d_acorr
        P_acorr_res[i,j] = P_acorr
        fwhm_acorr_res[i,j] = fwhm_acorr
    
# mean across iterations
k_fft_mean = np.nanmean(k_fft_res,0)
P_fft_mean = np.nanmean(P_fft_res,0)
d_acorr_mean = np.nanmean(d_acorr_res,0)
P_acorr_mean = np.nanmean(P_acorr_res,0)
fwhm_acorr_mean = np.nanmean(fwhm_acorr_res,0)
    
# std across iterations
k_fft_std = np.nanstd(k_fft_res,0)
P_fft_std = np.nanstd(P_fft_res,0)
d_acorr_std = np.nanstd(d_acorr_res,0)
P_acorr_std = np.nanstd(P_acorr_res,0)
fwhm_acorr_std = np.nanstd(fwhm_acorr_res,0)

"""
Example plots
"""
# get odc pattern
_, neural, _, _, _, _ = odc_2d(Nx_sim, Ny_sim, FOVx, FOVy, Nx_mri, Ny_mri, rho, delta, epsilon, 
                               theta, alpha, beta, fwhm_bold, fwhm_noise, a, b, theta, False)

# get pca
x_major, y_major, x_minor, y_minor = get_pca(neural)

# analyze fft
_, _, x1, y1 = analyze_fft(neural, FOVx, FOVy, x_major, y_major, 10, None)
_, _, x2, y2 = analyze_fft(neural, FOVx, FOVy, x_minor, y_minor, 10, None)

# fft
fig, ax = plt.subplots()
ax.plot(x1[:50], y1[:50])
ax.plot(x2[:50], y2[:50])
ax.set_xlabel("Spatial frequency in cycles/mm")
ax.set_ylabel("Power spectrum in a.u.")
ax.legend(["major axis","minor axis"])
fig.savefig(os.path.join(path_output,"fft_example.png"), bbox_inches="tight")
plt.show()

# analyze autocorrelation
_, _, _, x1, y1 = analyze_acorr(neural, FOVx, FOVy, x_major, y_major, 0.01, None)
_, _, _, x2, y2 = analyze_acorr(neural, FOVx, FOVy, x_minor, y_minor, 0.01, None)

# nac
fig, ax = plt.subplots()
ax.plot(x1, y1)
ax.plot(x2, y2)
ax.set_xlabel("Lag in mm")
ax.set_ylabel("NAC")
ax.legend(["major axis","minor axis"])
fig.savefig(os.path.join(path_output,"acorr_example.png"), bbox_inches="tight")
plt.show()

"""
Line plots
"""
# k_fft plot
fig, ax = plt.subplots()
ax.plot(phi, k_fft_mean, "r")
ax.fill_between(phi, k_fft_mean-k_fft_std, k_fft_mean+k_fft_std, alpha=0.2, facecolor='b')
ax.set_xlabel("Rotation angle in deg")
ax.set_ylabel("Peak spatial frequency in cycles/mm")
fig.savefig(os.path.join(path_output,"frequency_line.png"), bbox_inches="tight")
plt.show()

# P_fft plot
fig, ax = plt.subplots()
ax.plot(phi, P_fft_mean, "r")
ax.fill_between(phi, P_fft_mean-P_fft_std, P_fft_mean+P_fft_std, alpha=0.2, facecolor='b')
ax.set_xlabel("Rotation angle in deg")
ax.set_ylabel("Peak spatial frequency power in a.u.")
fig.savefig(os.path.join(path_output,"frequency_power_line.png"), bbox_inches="tight")
plt.show()

# d_acorr plot
fig, ax = plt.subplots()
ax.plot(phi, d_acorr_mean, "r")
ax.fill_between(phi, d_acorr_mean-d_acorr_std, d_acorr_mean+d_acorr_std, alpha=0.2, facecolor='b')
ax.set_xlabel("Rotation angle in deg")
ax.set_ylabel("Nearest neighbor distance in mm")
fig.savefig(os.path.join(path_output,"neighbor_distance_line.png"), bbox_inches="tight")
plt.show()

# P_acorr plot
fig, ax = plt.subplots()
ax.plot(phi, P_acorr_mean, "r")
ax.fill_between(phi, P_acorr_mean-P_acorr_std, P_acorr_mean+P_acorr_std, alpha=0.2, facecolor='b')
ax.set_xlabel("Rotation angle in deg")
ax.set_ylabel("Nearest neighbor power in a.u.")
fig.savefig(os.path.join(path_output,"neighbor_power_line.png"), bbox_inches="tight")
plt.show()

# fwhm plot
fig, ax = plt.subplots()
ax.plot(phi, fwhm_acorr_mean, "r")
ax.fill_between(phi, fwhm_acorr_mean-fwhm_acorr_std, fwhm_acorr_mean+fwhm_acorr_std, alpha=0.2, 
                facecolor='b')
ax.set_xlabel("Rotation angle in deg")
ax.set_ylabel("Column width in mm")
fig.savefig(os.path.join(path_output,"fwhm_line.png"), bbox_inches="tight")
plt.show()
  
"""
Polar plots
"""
phi_rad = phi.copy()
phi_rad = phi / 180 * np.pi
phi_rad = np.append(phi_rad,0) 

# fwhm
r = fwhm_acorr_mean.copy()
r = np.append(r, r[0])
r_err = fwhm_acorr_std.copy()
r_err = np.append(r_err, r_err[0])

fig = plt.figure()
ax = plt.axes(polar=True)
ax.plot(phi_rad, r, 'ro')
ax.errorbar(phi_rad, r, yerr=r_err, xerr=0, capsize=0, ecolor="g")
ax.set_xlabel("Major axis")
ax.set_ylabel("Minor axis", labelpad=30)
ax.set_rticks([1.0, 1.5, 2.0, 2.5])
ax.set_rlabel_position(-22.5)
ax.grid(True)
fig.savefig(os.path.join(path_output,"fwhm_polar.png"), bbox_inches="tight")
plt.show()

# neighbor distance
r = d_acorr_mean.copy()
r = np.append(r, r[0])
r_err = d_acorr_std.copy()
r_err = np.append(r_err, r_err[0])

fig = plt.figure()
ax = plt.axes(polar=True)
ax.plot(phi_rad, r, 'ro')
ax.errorbar(phi_rad, r, yerr=r_err, xerr=0, capsize=0, ecolor="g")
ax.set_xlabel("Major axis")
ax.set_ylabel("Minor axis", labelpad=30)
ax.set_rticks([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
ax.set_rlabel_position(-22.5)  
ax.grid(True)
fig.savefig(os.path.join(path_output,"neighbor_distance_polar.png"), bbox_inches="tight")
plt.show()

# neighbor power
r = P_acorr_mean.copy()
r = np.append(r, r[0])
r_err = P_acorr_std.copy()
r_err = np.append(r_err, r_err[0])

fig = plt.figure()
ax = plt.axes(polar=True)
ax.plot(phi_rad, r, 'ro')
ax.errorbar(phi_rad, r, yerr=r_err, xerr=0, capsize=0, ecolor="g")
ax.set_xlabel("Major axis")
ax.set_ylabel("Minor axis", labelpad=30)
ax.set_rticks([0.1, 0.2, 0.3, 0.4])
ax.set_rlabel_position(-22.5) 
ax.grid(True)
fig.savefig(os.path.join(path_output,"neighbor_power_polar.png"), bbox_inches="tight")
plt.show()

# spatial frequency
r = k_fft_mean.copy()
r = np.append(r, r[0])
r_err = k_fft_std.copy()
r_err = np.append(r_err, r_err[0])

fig = plt.figure()
ax = plt.axes(polar=True)
ax.plot(phi_rad, r, 'ro')
ax.errorbar(phi_rad, r, yerr=r_err, xerr=0, capsize=0, ecolor="g")
ax.set_xlabel("Major axis")
ax.set_ylabel("Minor axis", labelpad=30)
ax.set_rticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
ax.set_rlabel_position(-22.5)
ax.grid(True)
fig.savefig(os.path.join(path_output,"frequency_polar.png"), bbox_inches="tight")
plt.show()

# spatial frequency power
r = P_fft_mean.copy()
r = np.append(r, r[0])
r_err = P_fft_std.copy()
r_err = np.append(r_err, r_err[0])

fig = plt.figure()
ax = plt.axes(polar=True)
ax.plot(phi_rad, r, 'ro')
ax.errorbar(phi_rad, r, yerr=r_err, xerr=0, capsize=0, ecolor="g")
ax.set_xlabel("Major axis")
ax.set_ylabel("Minor axis", labelpad=30)
ax.set_rticks([1000, 2000, 3000, 4000])
ax.set_rlabel_position(-22.5)
ax.grid(True)
fig.savefig(os.path.join(path_output,"frequency_power_polar.png"), bbox_inches="tight")
plt.show()

"""
Logfile
"""
fileID = open(os.path.join(path_output,"result.txt"),"w")
fileID.write("fwhm_acorr_mean (minor): "+str(fwhm_acorr_mean[0])+"\n")
fileID.write("fwhm_acorr_std (minor): "+str(fwhm_acorr_std[0])+"\n")
fileID.write("fwhm_acorr_mean (major): "+str(fwhm_acorr_mean[int(np.argwhere(phi==90))])+"\n")
fileID.write("fwhm_acorr_std (major): "+str(fwhm_acorr_std[int(np.argwhere(phi==90))])+"\n")
fileID.write("d_acorr_mean (minor): "+str(d_acorr_mean[0])+"\n")
fileID.write("d_acorr_std (minor): "+str(d_acorr_std[0])+"\n")
fileID.write("d_acorr_mean (major): "+str(d_acorr_mean[int(np.argwhere(phi==90))])+"\n")
fileID.write("d_acorr_std (major): "+str(d_acorr_std[int(np.argwhere(phi==90))])+"\n")
fileID.write("P_acorr_mean (minor): "+str(P_acorr_mean[0])+"\n")
fileID.write("P_acorr_std (minor): "+str(P_acorr_std[0])+"\n")
fileID.write("P_acorr_mean (major): "+str(P_acorr_mean[int(np.argwhere(phi==90))])+"\n")
fileID.write("P_acorr_std (major): "+str(P_acorr_std[int(np.argwhere(phi==90))])+"\n")
fileID.write("k_fft_mean (minor): "+str(k_fft_mean[0])+"\n")
fileID.write("k_fft_std (minor): "+str(k_fft_std[0])+"\n")
fileID.write("k_fft_mean (major): "+str(k_fft_mean[int(np.argwhere(phi==90))])+"\n")
fileID.write("k_fft_std (major): "+str(k_fft_std[int(np.argwhere(phi==90))])+"\n")
fileID.write("P_fft_mean (minor): "+str(P_fft_mean[0])+"\n")
fileID.write("P_fft_std (minor): "+str(P_fft_std[0])+"\n")
fileID.write("P_fft_mean (major): "+str(P_fft_mean[int(np.argwhere(phi==90))])+"\n")
fileID.write("P_fft_std (major): "+str(P_fft_std[int(np.argwhere(phi==90))])+"\n")
fileID.close()
