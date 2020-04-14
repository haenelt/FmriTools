"""
Plots for autocorrelation analysis (experiment)

The purpose of the following script is to generate exemplary plots for the autocorrelation analysis
of ODC data. Data are averaged across hemispheres and two separate sessions. The following plots are 
produced:
    1. cortical depth vs. fwhm (phi=0; minor axis)
    2. cortical depth vs. k_fft (phi=0; minor axis)
    3. cortical depth vs. P_fft (phi=0; minor axis)
    4. phi vs. fwhm (averaged across cortical depth)
    5. phi vs. k_fft (averaged across cortical depth)
    6. phi vs. P_fft (averaged across cortical depth)
    7. fft for single layers (for each session)
    8. autocorrelation for single layers (for each session)

created by Daniel Haenelt
Date created: 18-04-2019
Last modified: 18-04-2019
"""
import os
import copy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib import rc

# load data from single sessions
lh1 = np.load("/home/daniel/mpi/conference/ohbm/ohbm_2019/poster/odc_exp/results/lh.GE_EPI4.npz")
lh2 = np.load("/home/daniel/mpi/conference/ohbm/ohbm_2019/poster/odc_exp/results/lh.GE_EPI5.npz")
rh1 = np.load("/home/daniel/mpi/conference/ohbm/ohbm_2019/poster/odc_exp/results/rh.GE_EPI4.npz")
rh2 = np.load("/home/daniel/mpi/conference/ohbm/ohbm_2019/poster/odc_exp/results/rh.GE_EPI5.npz")

path_output = "/home/daniel/mpi/conference/ohbm/ohbm_2019/poster/odc_exp/img_1p0"

""" do not edit below """

# make output folder
if not os.path.exists(path_output):
    os.mkdir(path_output)

# font parameters for plots
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

"""
cortical depth vs. fwhm, k_fft, P_fft along minor axis
"""

# fwhm
y1 = lh1["fwhm_acorr_phi"]
y2 = lh2["fwhm_acorr_phi"]
y3 = rh1["fwhm_acorr_phi"]
y4 = rh2["fwhm_acorr_phi"]

x_mean = np.nanmean([y1[:,0], y2[:,0], y3[:,0], y4[:,0]],0)
x_std = np.nanstd([y1[:,0], y2[:,0], y3[:,0], y4[:,0]],0)
y = np.linspace(0,1,10)

fig, ax = plt.subplots()
ax.plot(x_mean, y, color="red")
ax.fill_betweenx(y, x_mean - x_std, x_mean + x_std, alpha=0.2)
ax.set_xlabel("Column width in mm")
ax.set_ylabel("Cortical depth")
ax.set_title("ODCs along minor axis")
plt.yticks(y, np.round(y[::-1],1))
fig.savefig(os.path.join(path_output,"exp_fwhm_vs_depth.svg"), format='svg', bbox_inches='tight')
plt.show()

# k_fft
y1 = lh1["k_fft_phi"]
y2 = lh2["k_fft_phi"]
y3 = rh1["k_fft_phi"]
y4 = rh2["k_fft_phi"]

x_mean = np.nanmean([y1[:,0], y2[:,0], y3[:,0], y4[:,0]],0)
x_std = np.nanstd([y1[:,0], y2[:,0], y3[:,0], y4[:,0]],0)
y = np.linspace(0,1,10)

fig, ax = plt.subplots()
ax.plot(x_mean, y, color="red")
ax.fill_betweenx(y, x_mean - x_std, x_mean + x_std, alpha=0.2)
ax.set_xlabel("Peak spatial frequency in cycles/mm")
ax.set_ylabel("Cortical depth")
ax.set_title("ODCs along minor axis")
plt.yticks(y, np.round(y[::-1],1))
fig.savefig(os.path.join(path_output,"exp_column_kfft_vs_depth.svg"), format='svg', bbox_inches='tight')
plt.show()

# P_fft
y1 = lh1["P_fft_phi"]
y2 = lh2["P_fft_phi"]
y3 = rh1["P_fft_phi"]
y4 = rh2["P_fft_phi"]

x_mean = np.nanmean([y1[:,0], y2[:,0], y3[:,0], y4[:,0]],0)
x_std = np.nanstd([y1[:,0], y2[:,0], y3[:,0], y4[:,0]],0)
y = np.linspace(0,1,10)

fig, ax = plt.subplots()
ax.plot(x_mean, y, color="red")
ax.fill_betweenx(y, x_mean - x_std, x_mean + x_std, alpha=0.2)
ax.set_xlabel("Peak spatial frequency magnitude in a.u.")
ax.set_ylabel("Cortical depth")
ax.set_title("ODCs along minor axis")
plt.yticks(y, np.round(y[::-1],1))
fig.savefig(os.path.join(path_output,"exp_column_Pfft_vs_depth.svg"), format='svg', bbox_inches='tight')
plt.show()

"""
phi vs. fwhm, k_fft, P_fft
"""

# fwhm
y1 = lh1["fwhm_acorr_phi"]
y2 = lh2["fwhm_acorr_phi"]
y3 = rh1["fwhm_acorr_phi"]
y4 = rh2["fwhm_acorr_phi"]

phi = 10*np.arange(36)
y = np.concatenate((y1,y2,y3,y4),0)
y_mean = np.nanmean(y,0)
y_std = np.nanstd(y,0)

fig, ax = plt.subplots()
ax.plot(phi, y_mean, "r")
ax.fill_between(phi, y_mean-y_std, y_mean+y_std, alpha=0.2, facecolor='b')
ax.set_xlabel("Rotation angle in deg")
ax.set_ylabel("Column width in mm")
ax.set_title("ODCs along different projections")
fig.savefig(os.path.join(path_output,"exp_fwhm_vs_phi.svg"), format='svg', bbox_inches='tight')
plt.show()

# k_fft
y1 = lh1["k_fft_phi"]
y2 = lh2["k_fft_phi"]
y3 = rh1["k_fft_phi"]
y4 = rh2["k_fft_phi"]

phi = 10*np.arange(36)
y = np.concatenate((y1,y2,y3,y4),0)
y_mean = np.nanmean(y,0)
y_std = np.nanstd(y,0)

fig, ax = plt.subplots()
ax.plot(phi, y_mean, "r")
ax.fill_between(phi, y_mean-y_std, y_mean+y_std, alpha=0.2, facecolor='b')
ax.set_xlabel("Rotation angle in deg")
ax.set_ylabel("Peak spatial frequency in cycles/mm")
ax.set_title("ODCs along different projections")
fig.savefig(os.path.join(path_output,"exp_kfft_vs_phi.svg"), format='svg', bbox_inches='tight')
plt.show()

# Pfft
y1 = lh1["P_fft_phi"]
y2 = lh2["P_fft_phi"]
y3 = rh1["P_fft_phi"]
y4 = rh2["P_fft_phi"]

phi = 10*np.arange(36)
y = np.concatenate((y1,y2,y3,y4),0)
y_mean = np.nanmean(y,0)
y_std = np.nanstd(y,0)

fig, ax = plt.subplots()
ax.plot(phi, y_mean, "r")
ax.fill_between(phi, y_mean-y_std, y_mean+y_std, alpha=0.2, facecolor='b')
ax.set_xlabel("Rotation angle in deg")
ax.set_ylabel("Peak spatial frequency magnitude in a.u.")
ax.set_title("ODCs along different projections")
fig.savefig(os.path.join(path_output,"exp_Pfft_vs_phi.svg"), format='svg', bbox_inches='tight')
plt.show()

"""
FFT
"""
FOVx_mri = 148 # in mm
FOVy_mri = 148
Nx_mri = 148
Ny_mri = 148

# nyquist frequency
f_ny = np.sqrt(Nx_mri**2+Ny_mri**2) / (2 * np.sqrt(FOVx_mri**2+FOVy_mri**2))

# labels
labels = ["GM/WM"," "," "," "," "," "," "," "," ","GM/CSF"] # labels

for i in range(4):
    fig, ax = plt.subplots()
    axins = zoomed_inset_axes(ax, 2, loc=1, bbox_to_anchor= (10, 300)) # zoom-factor: 2.5, location: upper-left
    for j in range(10):
        if i == 0:
            y = lh1["y_fft_0"][j]
            x = lh1["x_fft_0"][j]
            name_label = "(LH; session 1)"
            plot_suffix = "lh_sess1"
        elif i == 1:
            y = lh2["y_fft_0"][j]
            x = lh2["x_fft_0"][j]
            name_label = "(LH; session 2)"
            plot_suffix = "lh_sess2"
        elif i == 2:
            y = rh1["y_fft_0"][j]
            x = rh1["x_fft_0"][j]
            name_label = "(RH; session 1)"
            plot_suffix = "rh_sess1"
        else:
            y = rh2["y_fft_0"][j]
            x = rh2["x_fft_0"][j]
            name_label = "(RH; session 2)"
            plot_suffix = "rh_sess2"
        y = y[x < 2]
        x = x[x < 2]
        ax.plot(x,y,label=labels[j],linewidth=1.0)
        axins.plot(x, y)
    ax.set_xlabel("Spatial frequency in cycles/mm")
    ax.set_ylabel("Normalized magnitude in a.u.")
    ax.set_title("Fourier spectrum along minor axis "+name_label, y=1.08)
    handles, labels = ax.get_legend_handles_labels()
    handles = [copy.copy(ha) for ha in handles ]
    [ha.set_linewidth(7) for ha in handles ] # set the linewidths to the copies
    ax.legend(handles[::-1], labels[::-1], loc=1, frameon=False, labelspacing=-0.3)
    ax.axvline(x=f_ny, ymin=0, ymax=1, linestyle='--', color='red')
    ax.text(f_ny-0.1,108,r"$f_{\mathrm{Nyquist}}$")
    x1, x2, y1, y2 = 0.1, 0.42, 0, 30 # specify the limits
    axins.set_xlim(x1, x2) # apply the x-limits
    axins.set_ylim(y1, y2) # apply the y-limits
    plt.yticks(visible=False)
    plt.xticks(visible=False)
    mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")
    fig.savefig(os.path.join(path_output,"exp_fft_example_"+plot_suffix+".svg"), format='svg', bbox_extra_artists=(axins,), bbox_inches='tight')
    plt.show()

"""
autocorrelation
"""

# labels
labels = ["GM/WM"," "," "," "," "," "," "," "," ","GM/CSF"]

for i in range(4):
    fig, ax = plt.subplots()
    axins = zoomed_inset_axes(ax, 2, loc=1, bbox_to_anchor= (10, 300)) # zoom-factor: 2.5, location: upper-left
    for j in range(10):
        if i == 0:
            y = lh1["y_acorr_0"][j]
            x = lh1["x_acorr_0"][j]
            name_label = "(LH; session 1)"
            plot_suffix = "lh_sess1"
        elif i == 1:
            y = lh2["y_acorr_0"][j]
            x = lh2["x_acorr_0"][j]
            name_label = "(LH; session 2)"
            plot_suffix = "lh_sess2"
        elif i == 2:
            y = rh1["y_acorr_0"][j]
            x = rh1["x_acorr_0"][j]
            name_label = "(RH; session 1)"
            plot_suffix = "rh_sess1"
        else:
            y = rh2["y_acorr_0"][j]
            x = rh2["x_acorr_0"][j]
            name_label = "(RH; session 2)"
            plot_suffix = "rh_sess2"
        ax.plot(x,y,label=labels[j],linewidth=1.0)
        axins.plot(x, y)
    ax.set_xlabel("Lag in mm")
    ax.set_ylabel("NAC in a.u.")
    ax.set_title("Autocorrelation along minor axis "+name_label)
    handles, labels = ax.get_legend_handles_labels()
    handles = [copy.copy(ha) for ha in handles ]
    [ha.set_linewidth(7) for ha in handles ] # set the linewidths to the copies
    ax.legend(handles[::-1], labels[::-1], loc=1, frameon=False, labelspacing=-0.3)
    x1, x2, y1, y2 = 1, 10, 0, 0.4 # specify the limits
    axins.set_xlim(x1, x2) # apply the x-limits
    axins.set_ylim(y1, y2) # apply the y-limits
    plt.yticks(visible=False)
    plt.xticks(visible=False)
    mark_inset(ax, axins, loc1=1, loc2=3, fc="none", ec="0.5")
    fig.savefig(os.path.join(path_output,"exp_acorr_example_"+plot_suffix+".svg"), format='svg', bbox_extra_artists=(axins,), bbox_inches='tight')
    plt.show()
