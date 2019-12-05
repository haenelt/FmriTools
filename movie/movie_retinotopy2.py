"""
Make retinotopy movie

The purpose of the following script is to make a movie to illustrate directly the demeaned time 
series (travelling wave) on the cortical surface. Optionally, the data will be additionally sampled 
on a regular grid.

created by Daniel Haenelt
Date created: 11-11-2019
Last modified: 04-12-2019
"""
import os
import re
import glob
import shutil as sh
import numpy as np
import nibabel as nb
from nighres.registration import apply_coordinate_mappings
from lib.processing.demean_time_series import demean_time_series
from lib.mapping import map2surface, map2stack

input_series = "/nobackup/actinium2/haenelt/V2STRIPES/p6/psf/multipol_14/udata.nii"
input_deformation = "/nobackup/actinium2/haenelt/V2STRIPES/p6/deformation/multipol/epi2orig.nii.gz"
input_surf = ["/nobackup/actinium2/haenelt/V2STRIPES/p6/anatomy/layer/lh.layer0",
              "/nobackup/actinium2/haenelt/V2STRIPES/p6/anatomy/layer/lh.layer1",
              "/nobackup/actinium2/haenelt/V2STRIPES/p6/anatomy/layer/lh.layer2",
              "/nobackup/actinium2/haenelt/V2STRIPES/p6/anatomy/layer/lh.layer3",
              "/nobackup/actinium2/haenelt/V2STRIPES/p6/anatomy/layer/lh.layer4",
              "/nobackup/actinium2/haenelt/V2STRIPES/p6/anatomy/layer/lh.layer5",
              "/nobackup/actinium2/haenelt/V2STRIPES/p6/anatomy/layer/lh.layer6",
              "/nobackup/actinium2/haenelt/V2STRIPES/p6/anatomy/layer/lh.layer7",
              "/nobackup/actinium2/haenelt/V2STRIPES/p6/anatomy/layer/lh.layer8",
              "/nobackup/actinium2/haenelt/V2STRIPES/p6/anatomy/layer/lh.layer9",
              ]
input_grid = ["/nobackup/actinium2/haenelt/V2STRIPES/p6/anatomy/ortho/lh.occip4.patch.flat.cmap.nii"]
path_output = "/data/pt_01880/multipol14_grid"

# paramteers
TR = 3
cutoff_highpass = 144
retinotopy_period = 48
start_vol = 11
end_vol = 4
interpolation = "linear" # can be linear or nearest
sigma_grid = 0.5

# path to SPM12 folder
pathSPM = "/data/pt_01880/source/spm12"
pathLIB = "/home/raid2/haenelt/projects/scripts/lib/preprocessing"

""" do not edit below """

# prefix
tmp = np.random.randint(0, 10, 5)
prefix = ''.join(str(i) for i in tmp)+'_'

# change to lib folder
os.chdir(pathLIB)

# make output folders
path_native = os.path.join(path_output,"native")
path_def = os.path.join(path_output,"def")
path_surf = os.path.join(path_output,"surf")
path_grid = os.path.join(path_output,"grid")

if not os.path.exists(path_output):
    os.mkdir(path_output)

if not os.path.exists(path_native):
    os.mkdir(path_native)
    
if not os.path.exists(path_def):
    os.mkdir(path_def)

if not os.path.exists(path_surf):
    os.mkdir(path_surf)

if not os.path.exists(path_grid):
    os.mkdir(path_grid)

# look for baseline corrected time series
os.system("matlab" + \
          " -nodisplay -nodesktop -r " + \
          "\"baseline_correction(\'{0}\', {1}, {2}, \'{3}\', \'{4}\'); exit;\"". \
          format(input_series, TR, cutoff_highpass, pathSPM, prefix))

# move baseline corrected time series to output folder
sh.move(os.path.join(os.path.dirname(input_series),prefix+os.path.basename(input_series)),
        os.path.join(path_native,"temp.nii"))

# demean time series
data = demean_time_series(os.path.join(path_native,"temp.nii"), write_output=False)
data_array = data.get_fdata()
    
# discard volumes at the beginning and at the end
if start_vol != 0:
   data_array = data_array[:,:,:,start_vol:]
   
if end_vol != 0:
    data_array = data_array[:,:,:,:-end_vol]

# squeeze retinotopy data to one period
period_length = int(retinotopy_period / TR)
retinotopy_array = np.zeros(np.append(data.header["dim"][1:4], period_length))
retinotopy_sort = np.mod(np.arange(0,np.shape(data_array)[3]), period_length)
for i in range(len(retinotopy_sort)):
    retinotopy_array[:,:,:,retinotopy_sort[i]] += data_array[:,:,:,i] 

# write single volumes
data.header["dim"][4] = 1
data.header["dim"][0] = 3
for i in range(np.shape(retinotopy_array)[3]):
    output = nb.Nifti1Image(retinotopy_array[:,:,:,i], data.affine, data.header)
    nb.save(output, os.path.join(path_native,str(i+1)+".nii")) 

# apply deformation    
for i in range(np.shape(retinotopy_array)[3]):
    apply_coordinate_mappings(os.path.join(path_native,str(i+1)+".nii"), 
                              input_deformation, 
                              interpolation=interpolation, 
                              padding='closest', 
                              save_data=True, 
                              overwrite=True, 
                              output_dir=path_def,
                              file_name=None,
                              )

# map to ana
for i in range(np.shape(retinotopy_array)[3]):
    for j in range(len(input_surf)):

        # hemisphere
        hemi = os.path.splitext(os.path.basename(input_surf[j]))[0]
        
        # sample on surface
        map2surface(input_surf[j], 
                    os.path.join(path_def,str(i+1)+"_def-img.nii.gz"), 
                    hemi, 
                    path_surf, 
                    input_white=None,
                    input_ind=None,
                    cleanup=True)

# sample on grid (optional)
if len(input_grid) > 0:
    for i in range(len(input_grid)):
        
        # hemisphere
        hemi = os.path.basename(input_grid[i])[:2]
        
        # sample on grid
        files = glob.glob(os.path.join(path_surf,'*'))
        for j in range(np.shape(retinotopy_array)[3]):
            stack = []
            nstack = []
            for k in range(len(files)):
                if re.findall(r''+hemi+'.'+str(j+1)+'_', files[k]): # quite nasty
                    stack.append(files[k])
                    nstack.append(re.findall(r'\d+', files[k])[-1])
                
            # string to number
            nstack = [int(x) for x in nstack]
            
            # get sorted list
            stack = [x for _,x in sorted(zip(nstack,stack), reverse=True)]
            
            # map to stack
            map2stack(stack, input_grid[i], sigma_grid, path_grid)

# delete temp
os.remove(os.path.join(path_native,"temp.nii"))