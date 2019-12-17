def calculate_equivolumetric_epi():
    """
    This function blas...
    Inputs:
        *input_white: bla
        *input_pial:
        *input_vol:
        *hemi:
        *
    
    created by Daniel Haenelt
    Date created: 17-12-2019
    Last modified: 17-12-2019
    """
#import os
import os
import subprocess
import numpy as np
import nibabel as nb
from nibabel.affines import apply_affine
from nibabel.freesurfer.io import read_geometry
from skimage import measure
from nighres.surface import probability_to_levelset
from nighres.laminar import volumetric_layering
from lib.utils.upsample_volume import upsample_volume
from lib.surface.upsample_surf_mesh import upsample_surf_mesh

input_white = "/home/raid2/haenelt/projects/GBB/test_data/lh.layer10_def"
input_pial = "/home/raid2/haenelt/projects/GBB/test_data/lh.layer0_def"
input_vol = "/home/raid2/haenelt/projects/GBB/test_data/mean_data.nii"
path_output = "/data/pt_01880/bla9"
hemi = "lh"

# parameters
dx = 0.4
dy = 0.4
dz = 0.4
niter = 2
n_start = 10
n_end = 5
n_layers = 2

""" do not edit below """

# make output folder
if not os.path.exists(path_output):
    os.makedirs(path_output)

# new filenames in output folder
res_white = os.path.join(path_output,hemi+".white")
res_pial = os.path.join(path_output,hemi+".pial")
res_vol = os.path.join(path_output,"epi_upsampled.nii")

# upsample volume
upsample_volume(input_vol, 
                res_vol,
                dxyz=[dx, dy, dz], 
                rmode="Cu")

# upsamples surface mesh
upsample_surf_mesh(input_white, 
                   res_white,
                   niter,
                   "linear")

upsample_surf_mesh(input_pial, 
                   res_pial,
                   niter, 
                   "linear")

# get affine ras2vox-tkr transformation to reference volume
transformation = subprocess.check_output(['mri_info', res_vol, '--{}'.format("ras2vox-tkr")]).decode()
num_transformation = [[float(x) for x in line.split()] for line in transformation.split('\n') if len(line)>0]
ras2vox_tkr = np.array(num_transformation)

# load surface
vtx_white, fac_white = read_geometry(res_white) 
vtx_pial, _ = read_geometry(res_pial)

# load volume
vol = nb.load(res_vol)

# apply ras2vox to coords
vtx_white = apply_affine(ras2vox_tkr, vtx_white).astype(int)
vtx_pial = apply_affine(ras2vox_tkr, vtx_pial).astype(int)

# surfaces to lines in volume
white_array = np.zeros(vol.header["dim"][1:4])
white_array[vtx_white[:,0],vtx_white[:,1],vtx_white[:,2]] = 1
white = nb.Nifti1Image(white_array, vol.affine, vol.header)   

pial_array = np.zeros(vol.header["dim"][1:4])
pial_array[vtx_pial[:,0],vtx_pial[:,1],vtx_pial[:,2]] = 1
pial = nb.Nifti1Image(pial_array, vol.affine, vol.header)

# lines to levelset
white_level = probability_to_levelset(white)
white_level_array = white_level["result"].get_fdata()

pial_level = probability_to_levelset(pial)
pial_level_array = pial_level["result"].get_fdata()

"""
make wm
"""
white_label_array = np.zeros_like(white_level_array)
white_label_array[white_level_array > pial_level_array] = 1
white_label_array[white_label_array != 1] = 0
white_label_array -= 1
white_label_array = np.abs(white_label_array).astype(int)
white_label_array = white_label_array - white_array
white_label_array[white_label_array < 0] = 0
white_label_array[:,:,:n_start] = 0
white_label_array[:,:,-n_end:] = 0
white_label_array = measure.label(white_label_array, connectivity=1)
white_label_array[white_label_array == 1] = 0
white_label_array[white_label_array > 0] = 1
white_label = nb.Nifti1Image(white_label_array, vol.affine, vol.header)

"""
make csf
"""
pial_label_array = np.zeros_like(pial_level_array)
pial_label_array[pial_level_array < white_level_array] = 1
pial_label_array[pial_label_array != 1] = 0
pial_label_array = np.abs(pial_label_array).astype(int)
pial_label_array = pial_label_array - pial_array
pial_label_array[pial_label_array < 0] = 0
pial_label_array[:,:,:n_start] = 0
pial_label_array[:,:,-n_end:] = 0
pial_label_array = measure.label(pial_label_array, connectivity=1)
pial_label_array[pial_label_array != 1] = 0
pial_label_array -= 1
pial_label_array = np.abs(pial_label_array)
pial_label_array[pial_array == 1] = 0
pial_label_array[:,:,:n_start] = 0
pial_label_array[:,:,-n_end:] = 0
pial_label = nb.Nifti1Image(pial_label_array, vol.affine, vol.header)

"""
make gm
"""
ribbon_label_array = pial_label_array.copy()
ribbon_label_array -= 1
ribbon_label_array = np.abs(ribbon_label_array)
ribbon_label_array = ribbon_label_array + white_label_array
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
                    file_name="epi")

# write niftis
nb.save(white, os.path.join(path_output,"wm_line.nii"))
nb.save(pial, os.path.join(path_output,"csf_line.nii"))
nb.save(white_label,os.path.join(path_output,"wm_label.nii"))
nb.save(pial_label,os.path.join(path_output,"csf_label.nii"))
nb.save(ribbon_label, os.path.join(path_output,"gm_label.nii"))