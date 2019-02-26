"""
Reliability measure between sessions. Data within a defined label is compared vertex-wise using 
Pearson correlation. Only a fraction of vertices is used to account for spatial covariance of nearby 
vertices. The p-value is estimated using a permutation analysis. A scatter and a Bland-Altman plot 
is created.

created by Daniel Haenelt
Date created: 17-11-2018
Last modified: 19-11-2018
"""
import os
import numpy as np
import nibabel as nb
import random
import matplotlib.pyplot as plt
from nibabel.freesurfer.io import read_label
from scipy.stats import pearsonr

# input paths
input_label = '/home/daniel/Schreibtisch/home/label/somewhere.label'
input_sess1 = '/home/daniel/Schreibtisch/home/flat/spmT/lh.spmT_left_right_GE_EPI2_nomotion_def_layer2.mgh'
input_sess2 = '/home/daniel/Schreibtisch/home/flat/spmT/lh.spmT_left_right_GE_EPI4_nomotion_def_layer2.mgh'
path_output = '/home/daniel/Schreibtisch'

# parameters
frac = 0.1
niter = 10000

""" do not edit below """

# load data
label = read_label(input_label).tolist()
sess1 = nb.load(input_sess1).get_fdata()
sess2 = nb.load(input_sess2).get_fdata()

# get the amount of data points
#ndata = np.round(frac * len(label)).astype(int)
ndata = 87 # 10 percent of data points within the first roi 

# randomly select ndata points in sess1 and sess2
label_shuffled = random.sample(label, len(label))
label_shuffled = label_shuffled[0:ndata]

# get correlation between sessions
x = sess1[label_shuffled][:,0,0]
y = sess2[label_shuffled][:,0,0]
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

# scatter plot with regression line
plt.figure()
plt.scatter(x, y)
plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)),'r')
plt.title(os.path.basename(input_label))
plt.xlabel('t-score (session 1)')
plt.ylabel('t-score (session 2)')
plt.savefig(os.path.join(path_output,os.path.splitext((os.path.basename(input_label)))[0]+'.png'))

# compare correlation coefficient to change level (permutation)
null_dist = []
for i in range(niter):
    # permute data points in sess2
    label_shuffled = random.sample(label, len(label))
    label_shuffled = label_shuffled[0:ndata]
    
    y = sess2[label_shuffled][:,0,0]
    null_dist.append(pearsonr(x,y)[0])
    

p = np.array(null_dist) > r[0]
p = p[p == True]
p = len(p) / niter

# print results
print('ROI: '+os.path.basename(input_label))
print('Number of points: '+str(ndata))
print('Number of iterations: '+str(niter))
print('Correlation coefficient: '+str(r[0]))
print('p-value: '+str(p))
