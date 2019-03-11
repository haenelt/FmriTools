"""
Reliability measure between sessions. Data within a defined label is compared vertex-wise using 
Pearson correlation. Only a fraction of vertices is used to account for spatial covariance of nearby 
vertices. The p-value is estimated using a permutation analysis. A scatter and a Bland-Altman plot 
is created.

created by Daniel Haenelt
Date created: 17-11-2018
Last modified: 11-03-2019
"""
import os
import numpy as np
import nibabel as nb
import random
import matplotlib.pyplot as plt
from nibabel.freesurfer.io import read_label, read_morph_data
from scipy.stats import pearsonr

# input paths
input_label = '/data/pt_01880/test/lh.stripes_in.label'
input_sess1 = '/data/pt_01880/V2STRIPES/p6/resting_state/alff_biopac/surf/lh.falff_z_def-img.nii_layer5_def.mgh'
input_sess2 = '/data/pt_01880/V2STRIPES/p6/colour/results/spmT/surf/lh.spmT_colour_bw_GE_EPI1_def_layer5.mgh'
path_output = '/data/pt_01880'

# parameters
frac = 0.1
niter = 1000

""" do not edit below """

# load data
label = read_label(input_label).tolist()

# if input file extension is not *.mgh, interprete as morphological file
if os.path.splitext(os.path.basename(input_sess1))[1] == ".mgh":
    sess1 = np.squeeze(nb.load(input_sess1).get_fdata())
else:
    sess1 = np.squeeze(read_morph_data(input_sess1))

# if input file extension is not *.mgh, interprete as morphological file
if os.path.splitext(os.path.basename(input_sess2))[1] == ".mgh":
    sess2 = np.squeeze(nb.load(input_sess2).get_fdata())
else:
    sess2 = np.squeeze(read_morph_data(input_sess2))

# get the amount of data points
ndata = np.round(frac * len(label)).astype(int)

# randomly select ndata points in sess1 and sess2
label_shuffled = random.sample(label, len(label))
label_shuffled = label_shuffled[0:ndata]

# get correlation between sessions
x = sess1[label_shuffled]
y = sess2[label_shuffled]
r = pearsonr(x, y)

# Bland Altman plot
mean = np.mean([x, y], axis=0)
diff = x - y  # difference between data1 and data2
md = np.mean(diff)  # mean of the difference
sd = np.std(diff, axis=0)  # Standard deviation of the difference

plt.scatter(mean, diff)
plt.axhline(md,           color='gray', linestyle='--')
plt.axhline(md + 1.96*sd, color='gray', linestyle='--')
plt.axhline(md - 1.96*sd, color='gray', linestyle='--')
plt.title('Bland-Altman Plot: '+os.path.basename(input_label))
plt.xlabel("Average of 2 sessions")
plt.ylabel("Difference between 2 sessions")
plt.savefig(os.path.join(path_output,os.path.splitext((os.path.basename(input_label)))[0]+'_bland_altman.png'))

# compare correlation coefficient to change level (permutation)
null_dist = []
for i in range(niter):
    # permute data points in sess2
    label_shuffled = random.sample(label, len(label))
    label_shuffled = label_shuffled[0:ndata]
    
    y_shuffle = sess2[label_shuffled]
    null_dist.append(pearsonr(x,y_shuffle)[0])
    

p = np.array(null_dist) > r[0]
p = p[p == True]
p = len(p) / niter

# print results
print('ROI: '+os.path.basename(input_label))
print('Number of points: '+str(ndata))
print('Number of iterations: '+str(niter))
print('Correlation coefficient: '+str(r[0]))
print('p-value: '+str(p))

# scatter plot with regression line
plt.figure()
plt.scatter(x, y)
plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)),'r')
plt.title(os.path.basename(input_label))
plt.xlabel('t-score (session 1)')
plt.ylabel('t-score (session 2)')
plt.ylabel('curvature')
plt.figtext(0.27, 0.8, "r = "+str(np.round(r[0],3))+", p = "+str(np.round(p, 3)), horizontalalignment="center")
plt.savefig(os.path.join(path_output,os.path.splitext((os.path.basename(input_label)))[0]+'.png'))