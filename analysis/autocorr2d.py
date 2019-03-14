import os
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nb
from scipy.signal import correlate2d
from skimage.measure import find_contours
from numpy.fft import fft2, fftshift

"""
Daran muss noch gearbeitet werden. Die Idee ist eine räumliche Autokorrelation zu machen, durch eine
partielle PCA die Prinzipalkomponenten des Hauptpeaks zu bestimmen (mittels Berechnung der 
Kovarianzmatrix) und die FFT Komponente entlang der kleinere Hauptachse zu analysieren.
Es wird auch der Rotationswinkel zur Hauptachse ausgegeben. Damit kann der patch rotiert werden, um
die kolumnare Architektur im Querschnitt am besten zu zeigen.
"""

# input file
input = "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p6/grid/rh.spmT_left_right_GE_EPI4_def-img_layer9_def_sigma0.5_grid.nii"
path_output = "/home/raid2/haenelt/Desktop/ismrm_workshop_analysis/p6/grid"

# load input nifti
data_img = nb.load(input)
data_array = data_img.get_fdata()

# get autocorrelation and FFT for each slice
data_corr = np.zeros_like(data_array)
data_fft = np.zeros_like(data_array)
for i in range(np.size(data_array,2)):
    data1 = ( data_array[:,:,i] - np.mean(data_array[:,:,i]) ) / ( np.std(data_array[:,:,i]) *np.shape(data_array)[0]*np.shape(data_array)[1] ) 
    data2 = ( data_array[:,:,i] - np.mean(data_array[:,:,i]) ) / ( np.std(data_array[:,:,i]) ) 
    data_corr[:,:,i] = correlate2d(data1,data2,mode='same',boundary='fill',fillvalue=0)
    data_fft[:,:,i] = np.abs(fftshift(fft2(data_array[:,:,i])))
    
# write output
filenameOUT = os.path.join(path_output,os.path.splitext(os.path.basename(input))[0]+'_corr.nii')
data_img.header["dim"][0] = 2
data_img.header["dim"][3] = 1
output = nb.Nifti1Image(data_corr, data_img.affine, data_img.header)
nb.save(output,filenameOUT)

filenameOUT = os.path.join(path_output,os.path.splitext(os.path.basename(input))[0]+'_fft.nii')
data_img.header["dim"][0] = 2
data_img.header["dim"][3] = 1
output = nb.Nifti1Image(data_fft, data_img.affine, data_img.header)
nb.save(output,filenameOUT)

#%%

# get contour of central peak
layer = 4
central_peak = data_corr[:,:,layer].copy()
central_peak[central_peak < 0.1] = 0
central_peak[central_peak != 0] = 1
central_peak = find_contours(central_peak, level=0.1)

# define central peak contour as contour with longest outline (pretty ugly but is sufficient for the moment)
contour_length = 0
contour_ind = 0
for i in range(len(central_peak)):
    if len(central_peak[i]) > contour_length:
        contour_length = len(central_peak[i])
        contour_ind = i

# get x and y coordinates of the contour
x = []
y = []
for i in range(contour_length):
    x.append(central_peak[contour_ind][i][0].astype(int))
    y.append(central_peak[contour_ind][i][1].astype(int))

# define ellipse image from central peak contour
ellipse = np.zeros_like(data_corr[:,:,layer])
ellipse[x,y] = 1

# subtract mean from each dimension
x = np.array(x)
y = np.array(y)
x_mean = np.mean(x)
y_mean = np.mean(y)
x = ( x - x_mean )
y = ( y - y_mean )
coords = np.vstack([x, y])

# get covariance matrix
cov = np.cov(coords)
evals, evecs = np.linalg.eig(cov)

# sort eigenvalues in decreasing order
sort_indices = np.argsort(evals)[::-1]
x_v1, y_v1 = evecs[:, sort_indices[0]]  # Eigenvector with largest eigenvalue
x_v2, y_v2 = evecs[:, sort_indices[1]]

# plot the principal components
scale = 200
X_major = np.linspace(x_v1*-scale*2 + y_mean, x_v1*scale*2 + y_mean,1000)
Y_major = np.linspace(y_v1*-scale*2 + x_mean, y_v1*scale*2 + x_mean,1000)
X_minor = np.linspace(x_v2*-scale + y_mean, x_v2*scale + y_mean,1000)
Y_minor = np.linspace(y_v2*-scale + x_mean, y_v2*scale + x_mean,1000)

fig, ax = plt.subplots()
ax.imshow(data_corr[:,:,4])
ax.plot(X_major, Y_major, color='red')
ax.plot(X_minor, Y_minor, color='blue')
plt.show()

# get coordinates of minor axis as nearest neighbor of image pixels
line = []
for i in range(len(X_minor)):
    if X_minor[i] < 0:
        continue
    elif Y_minor[i] < 0:
        continue
    elif X_minor[i] > np.shape(data_fft[:,:,4])[0] or Y_minor[i] > np.shape(data_fft[:,:,4])[1]:
        continue
    else:
        x = np.round(X_minor[i]).astype(int)
        y = np.round(Y_minor[i]).astype(int)
        line.append(data_fft[y,x,4])

# plot line
line = np.array(line)
coord = np.linspace(0,len(line)*0.25,len(line))
plt.plot(coord,line)

# get the rotation angle
theta = np.tanh((x_v1)/(y_v1))/np.pi*180

# write output
filenameOUT = os.path.join("/nobackup/actinium1/haenelt",os.path.splitext(os.path.basename(input))[0]+'_corr.nii')
data_img.header["dim"][0] = 2
data_img.header["dim"][3] = 1
output = nb.Nifti1Image(data_corr, data_img.affine, data_img.header)
nb.save(output,filenameOUT)
