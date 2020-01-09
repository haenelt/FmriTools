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
input_white_surf = ["/home/daniel/Schreibtisch/temp/lh.layer10_def",
                    "/home/daniel/Schreibtisch/temp/rh.layer10_def",
                    ]

input_pial_surf = ["/home/daniel/Schreibtisch/temp/lh.layer0_def",
                   "/home/daniel/Schreibtisch/temp/rh.layer0_def",
                   ]

input_white_ind = ["/home/daniel/Schreibtisch/temp/lh.layer10_ind.txt",
                   "/home/daniel/Schreibtisch/temp/rh.layer10_ind.txt",
                   ]

input_pial_ind = ["/home/daniel/Schreibtisch/temp/lh.layer0_ind.txt",
                  "/home/daniel/Schreibtisch/temp/rh.layer0_ind.txt",
                  ]

# output directory
path_output = ""

""" do not edit below """

for i in range(len(input_white_surf)):
    
    # match vertex numbers
    match_vertex_number(input_white_surf[i], 
                        input_pial_surf[i], 
                        input_white_ind[i], 
                        input_pial_ind[i],
                        path_output)