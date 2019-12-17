# all renames in output folder
# make print commands
# clean folder
import os
import subprocess
import numpy as np
import nibabel as nb
from nibabel.affines import apply_affine
from nibabel.freesurfer.io import read_geometry
from nighres.surface import probability_to_levelset
from lib.utils import upsample_volume
from lib.segmentation import upsample_surf_mesh
from skimage import measure
from nighres.laminar import volumetric_layering

input_white = "/home/daniel/projects/GBB/test_data/lh.layer10_def"
input_pial = "/home/daniel/projects/GBB/test_data/lh.layer0_def"
input_vol = "/home/daniel/projects/GBB/test_data/mean_data.nii"
path_output = "/home/daniel/Schreibtisch/test"
hemi = "lh"

# parameters
dx = 0.4
dy = 0.4
dz = 0.4
niter = 2
n_start = 10
n_end = 5
n_layers = 3

""" do not edit below """

# make output folder
if not os.path.exists(path_output):
    os.makedirs(path_output)

# upsample volume
upsample_volume(input_vol, 
                os.path.join(path_output,"vol_upsampled.nii"), 
                dxyz=[dx, dy, dz], 
                rmode="Cu")

# upsamples surface mesh
upsample_surf_mesh(input_white, 
                   os.path.join(path_output,hemi+".white"),
                   niter, 
                   "linear")

upsample_surf_mesh(input_pial, 
                   os.path.join(path_output,hemi+".pial"),
                   niter, 
                   "linear")

# new filenames in output folder
input_white = os.path.join(path_output,hemi+".white")
input_pial = os.path.join(path_output,hemi+".pial")
input_vol = os.path.join(path_output,"vol_upsampled.nii")

# get affine ras2vox-tkr transformation to reference volume
transformation = subprocess.check_output(['mri_info', input_vol, '--{}'.format("ras2vox-tkr")]).decode()
num_transformation = [[float(x) for x in line.split()] for line in transformation.split('\n') if len(line)>0]
ras2vox_tkr = np.array(num_transformation)

# load surface
vtx_white, fac_white = read_geometry(input_white) 
vtx_pial, _ = read_geometry(input_pial)

# load volume
vol = nb.load(input_vol)

# apply ras2vox to coords
vtx_white = apply_affine(ras2vox_tkr, vtx_white).astype(int)
vtx_pial = apply_affine(ras2vox_tkr, vtx_pial).astype(int)

# surfaces to volume
white_array = np.zeros(vol.header["dim"][1:4])
white_array[vtx_white[:,0],vtx_white[:,1],vtx_white[:,2]] = 1
pial_array = np.zeros(vol.header["dim"][1:4])
pial_array[vtx_pial[:,0],vtx_pial[:,1],vtx_pial[:,2]] = 1

white = nb.Nifti1Image(white_array, vol.affine, vol.header)   
pial = nb.Nifti1Image(pial_array, vol.affine, vol.header)

# lines to levelset
white_level = probability_to_levelset(white)
pial_level = probability_to_levelset(pial)

white_level_array = white_level["result"].get_fdata()
pial_level_array = pial_level["result"].get_fdata()

white_final_array = np.zeros_like(white_level_array)
white_final_array[white_level_array > pial_level_array] = 1
white_final_array[white_final_array != 1] = 0
white_final_array -= 1
white_final_array = np.abs(white_final_array).astype(int)

white_final = nb.Nifti1Image(white_final_array, vol.affine, vol.header)

white_final_array2 = white_final_array - white_array
white_final_array2[white_final_array2 < 0] = 0



white_final_array2[:,:,:n_start] = 0
white_label_array = measure.label(white_final_array2, connectivity=1)
white_label_array[white_label_array == 1] = 0
white_label_array[white_label_array > 0] = 1


white_label = nb.Nifti1Image(white_label_array, vol.affine, vol.header)

"""
make pial surface
"""
pial_final_array = np.zeros_like(pial_level_array)
pial_final_array[pial_level_array < white_level_array] = 1
pial_final_array[pial_final_array != 1] = 0
#pial_final_array -= 1
pial_final_array = np.abs(pial_final_array).astype(int)

pial_final = nb.Nifti1Image(pial_final_array, vol.affine, vol.header)

pial_final_array2 = pial_final_array - pial_array
pial_final_array2[pial_final_array2 < 0] = 0

pial_final_array2[:,:,:n_start] = 0
pial_final_array2[:,:,-n_end:] = 0

pial_label_array = measure.label(pial_final_array2, connectivity=1)
pial_label_array[pial_label_array != 1] = 0
pial_label = nb.Nifti1Image(pial_label_array, vol.affine, vol.header)

"""
make ribbon
"""
ribbon_label_array = pial_label_array + white_label_array
ribbon_label_array -= 1
ribbon_label_array = np.abs(ribbon_label_array).astype(int)

ribbon_label_array[:,:,:n_start] = 0
ribbon_label_array[:,:,-n_end:] = 0

ribbon_label = nb.Nifti1Image(ribbon_label_array, vol.affine, vol.header)

"""
layers
"""
csf_level = probability_to_levelset(pial_label)
wm_level = probability_to_levelset(white_label)

volumetric_layering(wm_level["result"], 
                    csf_level["result"], 
                    n_layers=n_layers, 
                    topology_lut_dir=None,
                    save_data=True, 
                    overwrite=True, 
                    output_dir=path_output, 
                    file_name="layer")

#%%

nb.save(white_final,os.path.join(path_output,"white_final.nii"))
nb.save(white, os.path.join(path_output,"white_line.nii"))
nb.save(pial, os.path.join(path_output,"pial_line.nii"))
nb.save(white_final,os.path.join(path_output,"white_final.nii"))
nb.save(white_label,os.path.join(path_output,"white_label.nii"))
nb.save(pial_label,os.path.join(path_output,"pial_label.nii"))
nb.save(ribbon_label, os.path.join(path_output,"ribbon_label.nii"))
