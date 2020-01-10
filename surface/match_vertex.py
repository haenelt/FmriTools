"""
Match vertex numbers

The purpose of the following script is to match the number of vertices between deformed white and
pial surfaces. Vertices are removed which are not found in both surfaces. After removing, faces 
vertex-in-orig indices are updated.

created by Daniel Haenelt
Date created: 10-12-2019
Last modified: 10-12-2019
"""
from lib.surface import match_vertex_number

# input files
input_white_surf = ["/data/pt_01880/odc_temp/surface/def/lh.layer10_def2_smooth",
                   "/data/pt_01880/odc_temp/surface/def/rh.layer10_def2_smooth",
                   ]

input_pial_surf = ["/data/pt_01880/odc_temp/surface/def/lh.layer0_def2_smooth",
                   "/data/pt_01880/odc_temp/surface/def/rh.layer0_def2_smooth",
                   ]

input_white_ind = ["/data/pt_01880/odc_temp/surface/def/lh.layer10_def2_ind",
                  "/data/pt_01880/odc_temp/surface/def/rh.layer10_def2_ind",
                  ]

input_pial_ind = ["/data/pt_01880/odc_temp/surface/def/lh.layer0_def2_ind",
                  "/data/pt_01880/odc_temp/surface/def/rh.layer0_def2_ind",
                  ]

# output directory
path_output = "/data/pt_01880/odc_temp/surface/match"

""" do not edit below """

for i in range(len(input_white_surf)):
    
    # match vertex numbers
    match_vertex_number(input_white_surf[i], 
                        input_pial_surf[i], 
                        input_white_ind[i], 
                        input_pial_ind[i],
                        path_output)