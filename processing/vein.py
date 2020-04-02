"""
Venous mask

A binary mask of venous voxels is created from a set of time series. First, a baseline correction of
each time series is applied if not done before (i.e., if no file with prefix b is found). Venous 
voxels are classified in the following way: Mean epi and mean tSNR across all time series are 
thresholded using arbitrary thresholds (<epi_threshold>, <tsnr_threshold>) considering only low mean
and low tsnr voxels. The thresholds have to be adjusted individually. The final mask is created by 
intersecting both masks. Both mean epi and mean tsnr are saved for single runs and the whole 
session. Additionally, a text file is written containing the used parameters. The method is taken 
from Kashyap et al., 2017.

created by Daniel Haenelt
Date created: 06-12-2018             
Last modified: 13-09-2019
"""
import os
import datetime
import numpy as np
import nibabel as nb

# input data
img_input = ["/data/pt_01880/temp_gbb/data/udata.nii",
             ]

# path to SPM12 folder
pathSPM = "/data/pt_01880/source/spm12"
pathLIB = "/data/hu_haenelt/projects/scripts/lib"

# parameters
TR = 3 # repetition time in s
cutoff_highpass = 180 # cutoff in s for baseline correction
epi_threshold = 1200
tsnr_threshold = 12

""" do not edit below """

# change to lib folder
os.chdir(os.path.join(pathLIB,"preprocessing"))

# prepare path and filename
path = []
file = []
for i in range(len(img_input)):
    path.append(os.path.split(img_input[i])[0])
    file.append(os.path.split(img_input[i])[1])

# output folder is taken from the first entry of the input list and set into 
if len(img_input) > 1:
    path_output = os.path.join(os.path.dirname(path[0]),"vein","native")
else:
    path_output = os.path.join(path[0],"vein","native")
if not os.path.exists(path_output):
    os.makedirs(path_output)

# output filenames
name_epi = os.path.join(path_output,"mean_epi.nii")
name_tsnr = os.path.join(path_output,"mean_tsnr.nii")

# get image header information
data_img = nb.load(img_input[0])
data_img.header["dim"][0] = 3
data_img.header["dim"][4] = 1
header = data_img.header
affine = data_img.affine

# get image dimension
dim = data_img.header["dim"][1:4]

# mean epi of each time series
for i in range(len(path)):
    
    # look for baseline corrected time series
    if not os.path.isfile(os.path.join(path[i],"b"+file[i])):
        os.system("matlab" + \
                  " -nodisplay -nodesktop -r " + \
                  "\"baseline_correction(\'{0}\', {1}, {2}, \'{3}\'); exit;\"". \
                  format(img_input[i], TR, cutoff_highpass, pathSPM))

    if not os.path.isfile(os.path.join(path[i],"mean_b"+file[i])):
        # open baseline corrected time series
        data_img = nb.load(os.path.join(path[i],"b"+file[i]))
        data_array = data_img.get_fdata()

        # calculate mean of time series
        data_array_mean = np.mean(data_array, axis=3)

        # write output
        output = nb.Nifti1Image(data_array_mean, affine, header)
        fileOUT = os.path.join(path[i],"mean_b"+file[i])
        nb.save(output,fileOUT)
    
    if not os.path.isfile(os.path.join(path[i],"tsnr_b"+file[i])):
        # open baseline corrected time series
        data_img = nb.load(os.path.join(path[i],"b"+file[i]))
        data_array = data_img.get_fdata()
    
        # calculate mean and std of time series
        data_array_mean = np.mean(data_array, axis=3)
        data_array_std = np.std(data_array, axis=3)
        data_array_std[data_array_std == 0] = np.nan
    
        # tsnr of time series
        tsnr_array = data_array_mean / data_array_std
        tsnr_array[np.isnan(tsnr_array)] = 0
    
        # write output
        output = nb.Nifti1Image(tsnr_array, affine, header)
        fileOUT = os.path.join(path[i],"tsnr_b"+file[i])
        nb.save(output,fileOUT)

# mean epi across time series
if not os.path.isfile(name_epi):
    mean_epi_array = np.zeros(dim)  
    for i in range(len(path)):
        # load mean
        epi_img = nb.load(os.path.join(path[i],"mean_b"+file[i]))   
        epi_array = epi_img.get_fdata()    

        # sum mean arrays
        mean_epi_array += epi_array

    # mean epi
    mean_epi_array = mean_epi_array / len(path)

    # write output
    output = nb.Nifti1Image(mean_epi_array, affine, header)
    nb.save(output,name_epi)
else:
    mean_epi_img = nb.load(name_epi)
    mean_epi_array = mean_epi_img.get_fdata()

# mean tsnr across time series   
if not os.path.isfile(name_tsnr):
    mean_tsnr_array = np.zeros(dim)
    for i in range(len(path)):
        # load tsnr
        tsnr_img = nb.load(os.path.join(path[i],"tsnr_b"+file[i]))
        tsnr_array = tsnr_img.get_fdata()
    
        # sum mean arrays
        mean_tsnr_array += tsnr_array

    # mean tsnr
    mean_tsnr_array = mean_tsnr_array / len(path)

    output = nb.Nifti1Image(mean_tsnr_array, affine, header)
    nb.save(output,name_tsnr)
else:
    mean_tsnr_img = nb.load(name_tsnr)
    mean_tsnr_array = mean_tsnr_img.get_fdata()
   
# define masks
epi_mask_array = mean_epi_array.copy()
epi_mask_array[epi_mask_array > epi_threshold] = 0
epi_mask_array[epi_mask_array != 0] = 1

tsnr_mask_array = mean_tsnr_array.copy()
tsnr_mask_array[tsnr_mask_array > tsnr_threshold] = 0
tsnr_mask_array[tsnr_mask_array != 0] = 1

# combine masks
mask_array = epi_mask_array + tsnr_mask_array
mask_array[mask_array != 2] = 0
mask_array[mask_array != 0] = 1

# write output
output = nb.Nifti1Image(epi_mask_array, affine, header)
fileOUT = os.path.join(path_output,"mask_epi.nii")
nb.save(output,fileOUT)

output = nb.Nifti1Image(tsnr_mask_array, affine, header)
fileOUT = os.path.join(path_output,"mask_tsnr.nii")
nb.save(output,fileOUT)

output = nb.Nifti1Image(mask_array, affine, header)
fileOUT = os.path.join(path_output,"vein.nii")
nb.save(output,fileOUT)

# write log
fileID = open(os.path.join(path_output,"vein_info.txt"),"a")
fileID.write("script executed: "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
fileID.write("cutoff_highpass: "+str(cutoff_highpass)+"\n")
fileID.write("mean epi threshold: "+str(epi_threshold)+"\n")
fileID.write("tsnr threshold: "+str(tsnr_threshold)+"\n")
fileID.close()
