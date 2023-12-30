# -*- coding: utf-8 -*-
"""
Map2ana

In the following script, source data is transformed using a deformation field
and the transformed data is mapped onto a surface mesh in the target space. The
deformation of the input data is optional. The script needs an installation of
freesurfer.

"""

import os

from nighres.registration import apply_coordinate_mappings

from ..io.filename import get_filename
from ..mapping.map2surface import map2surface

# input
file_in = [
    "/data/pt_01983/func/meridian/results/psc/native/psc_vm_on_hm_on_GE_EPI1.nii",
]

surf_in = [
    "/data/pt_01983/func/anatomy/layer_mpm/Session1/lh.layer10_mpm",
    "/data/pt_01983/func/anatomy/layer_mpm/Session1/rh.layer10_mpm",
]

# parameters
deformation_in1 = "/data/pt_01983/func/deformation/meridian/GE_EPI1/epi2ana.nii.gz"
deformation_in2 = "/data/pt_01983/func/deformation/MPM/sess1/mp2rage_2_mpm.nii.gz"
path_output = "/data/pt_01983/func/meridian/results/psc"
interpolation = "linear"  # can be linear or nearest

# do not edit below

# output folders
path_def = os.path.join(path_output, "def")
path_surf = os.path.join(path_output, "surf")

for i in range(len(file_in)):
    # apply deformation
    if deformation_in1 is not None:
        apply_coordinate_mappings(
            file_in[i],
            deformation_in1,
            deformation_in2,
            interpolation=interpolation,
            padding="closest",
            save_data=True,
            overwrite=True,
            output_dir=path_def,
            file_name=None,
        )

    # get filename of deformed file
    if deformation_in1 is not None:
        _, name_file, _ = get_filename(file_in[i])
        filename_def = os.path.join(path_def, name_file + "_def-img.nii")
    else:
        filename_def = file_in[i]

    # map to ana
    for j in range(len(surf_in)):
        # sample on surface
        map2surface(
            surf_in[j],
            filename_def,
            True,
            path_surf,
            interp_method="nearest",
            input_surf_target=None,
            input_ind=None,
            cleanup=True,
        )
