# -*- coding: utf-8 -*-
"""
Make layer

This script computes equidistant layers between a white and a pial surface. The
following steps are done:
    (1) match vertex numbers of both surfaces
    (2) extract main component to remove unconnected parts
    (3) inflate white surface for visualization
    (4) compute meshlines
    (5) create layers

"""

import os

import numpy as np
from nibabel.freesurfer.io import read_geometry, write_geometry

from ..layer.get_meshlines import get_meshlines
from ..surface.extract_main_component import extract_main_component
from ..surface.inflate_surf_mesh import inflate_surf_mesh
from ..surface.match_vertex_number import match_vertex_number

# input
mesh = dict()
mesh["white"] = dict()
mesh["pial"] = dict()

ind = dict()
ind["white_1"] = dict()  # index file from first deformation
ind["white_2"] = dict()  # index file from second deformation
ind["pial_1"] = dict()
ind["pial_2"] = dict()

mesh["white"]["lh"] = "/data/pt_01880/test/gbb/lh.white_def2_refined"
mesh["white"]["rh"] = "/data/pt_01880/test/gbb/rh.white_def2_refined"
mesh["pial"]["lh"] = "/data/pt_01880/test/gbb/lh.pial_def2_refined"
mesh["pial"]["rh"] = "/data/pt_01880/test/gbb/rh.pial_def2_refined"

ind["white_1"][
    "lh"
] = "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/dense_epi/lh.white_def2_ind"
ind["white_1"][
    "rh"
] = "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/dense_epi/rh.white_def2_ind"
ind["white_2"]["lh"] = "/data/pt_01880/test/gbb/lh.white_def2_refined_ind.txt"
ind["white_2"]["rh"] = "/data/pt_01880/test/gbb/rh.white_def2_refined_ind.txt"
ind["pial_1"][
    "lh"
] = "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/dense_epi/lh.pial_def2_ind"
ind["pial_1"][
    "rh"
] = "/data/pt_01880/Experiment2_Rivalry/p1/anatomy/dense_epi/rh.pial_def2_ind"
ind["pial_2"]["lh"] = "/data/pt_01880/test/gbb/lh.pial_def2_refined_ind.txt"
ind["pial_2"]["rh"] = "/data/pt_01880/test/gbb/rh.pial_def2_refined_ind.txt"

path_output = "/data/pt_01880/test/blablabla2"
niter_inflate = 30
nlayer = 11

# do not edit below

# set output basenames
basename_white = "white"
basename_pial = "pial"

# make folders
path_dense = os.path.join(path_output, "dense_refined")
path_layer = os.path.join(path_output, "layer")

if not os.path.exists(path_output):
    os.makedirs(path_output)

if not os.path.exists(path_dense):
    os.makedirs(path_dense)
else:
    raise FileExistsError("Dense folder already exists!")

if not os.path.exists(path_layer):
    os.makedirs(path_layer)
else:
    raise FileExistsError("Layer folder already exists!")

hemi = ["lh", "rh"]
for i in range(len(hemi)):
    # load surfaces
    vtx_white, fac = read_geometry(mesh["white"][hemi[i]])
    vtx_pial, _ = read_geometry(mesh["pial"][hemi[i]])

    # load index files
    ind_white = np.loadtxt(ind["white_1"][hemi[i]]).astype(int)
    ind_white2 = np.loadtxt(ind["white_2"][hemi[i]]).astype(int)
    ind_white = ind_white[ind_white2]

    ind_pial = np.loadtxt(ind["pial_1"][hemi[i]]).astype(int)
    ind_pial2 = np.loadtxt(ind["pial_2"][hemi[i]]).astype(int)
    ind_pial = ind_pial[ind_pial2]

    # match vertex numbers
    vtx_white, vtx_pial, fac, _ = match_vertex_number(
        vtx_white, vtx_pial, fac, ind_white, ind_pial
    )

    file_white = os.path.join(path_dense, hemi[i] + "." + basename_white + "_match")
    file_pial = os.path.join(path_dense, hemi[i] + "." + basename_pial + "_match")

    write_geometry(file_white, vtx_white, fac)
    write_geometry(file_pial, vtx_pial, fac)

    # extract main component
    extract_main_component(file_white, file_white + "_final")
    extract_main_component(file_pial, file_pial + "_final")

    # inflation
    inflate_surf_mesh(
        file_white + "_final",
        file_white + "_final_inflate" + str(niter_inflate),
        niter_inflate,
    )

    # meshlines
    vtx_white, fac = read_geometry(file_white + "_final")
    vtx_pial, _ = read_geometry(file_pial + "_final")

    vtx_line, fac_line = get_meshlines(vtx_pial, vtx_white)
    file_line = os.path.join(path_dense, hemi[i] + ".mesh_lines")
    write_geometry(file_line, vtx_line, fac_line)

    # layer
    for j in range(nlayer):
        vtx_layer = vtx_white + j / (nlayer - 1) * (vtx_pial - vtx_white)
        file_layer = os.path.join(path_layer, hemi[i] + ".layer_" + str(j))
        write_geometry(file_layer, vtx_layer, fac)
