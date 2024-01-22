# -*- coding: utf-8 -*-
"""Make layer

This script computes equidistant layers between a white and a pial surface. The
following steps are done:
    (1) match vertex numbers of both surfaces
    (2) extract main component to remove unconnected parts
    (3) inflate white surface for visualization
    (4) compute meshlines
    (5) create layers

"""

import os
from argparse import SUPPRESS, ArgumentParser

import numpy as np
from nibabel.freesurfer.io import read_geometry, write_geometry

from ..io.filename import get_filename
from ..segmentation.layer import get_meshlines
from ..segmentation.surf import (
    extract_main_component,
    inflate_surf_mesh,
    match_vertex_number,
)

# parameters
HEMI = ["lh", "rh"]


def _get_parser():
    """Parse command line inputs.

    Returns
    -------
    parser.parse_args() : argparse dict
    """
    # Disable default help
    parser = ArgumentParser(add_help=False)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    # Add back help
    optional.add_argument(
        "--help",
        action="help",
        default=SUPPRESS,
        help="Show this help message and exit.",
    )
    required.add_argument(
        "--white",
        dest="surf_white",
        type=str,
        help="File name of white surface.",
        required=True,
    )
    required.add_argument(
        "--pial",
        dest="surf_pial",
        type=str,
        help="File name of pial surface.",
        required=True,
    )
    required.add_argument(
        "--ind_white1",
        dest="ind_white1",
        type=str,
        help="File name of first white label.",
        required=True,
    )
    required.add_argument(
        "--ind_white2",
        dest="ind_white2",
        type=str,
        help="File name of second white label.",
        required=True,
    )
    required.add_argument(
        "--ind_pial1",
        dest="ind_pial1",
        type=str,
        help="File name of first pial label.",
        required=True,
    )
    required.add_argument(
        "--ind_pial2",
        dest="ind_pial2",
        type=str,
        help="File name of second pial label.",
        required=True,
    )
    required.add_argument(
        "--out",
        dest="path_output",
        type=str,
        help="Directory where output is written.",
        required=True,
    )
    optional.add_argument(
        "--inflate",
        dest="n_iter_inflate",
        type=int,
        help=("Number of inflation iterations (default: %(default)s)."),
        default=30,
    )
    optional.add_argument(
        "--layer",
        dest="n_layer",
        type=int,
        help=("Number of equivolumetric layers (default: %(default)s)."),
        default=11,
    )

    return parser


def layer_workflow(
    surf_white,
    surf_pial,
    ind_white1,
    ind_white2,
    ind_pial1,
    ind_pial2,
    path_output,
    n_iter_inflate=30,
    n_layer=11,
):
    """Workflow for layer generation.

    Parameters
    ----------
    surf_white : str
        File name of white surface.
    surf_pial : str
        File name of pial surface.
    ind_white1 : str
        Set of vertex indices from first deformation (registration).
    ind_white2 : str
        Set of vertex indices from second deformation (GBB).
    ind_pial1 : str
        Set of vertex indices from first deformation (registration).
    ind_pial2 : str
        Set of vertex indices from second deformation (GBB).
    path_output : str
        Directory where output is written.
    n_iter_inflate : int, optional
        Number of inflation iterations, by default 30.
    n_layer : int, optional
        Number of layers, by default 11.

    Raises
    ------
    ValueError
        If hemisphere prefix is not found.
    ValueError
        If hemispheres of white and pial surfaces do not match.
    FileExistsError
        If output folder dense already exists.
    FileExistsError
        If output folder layer already exists.
    """
    # get hemisphere from input file
    _, hemi, _ = get_filename(surf_white)
    _, hemi_p, _ = get_filename(surf_pial)
    if (hemi not in HEMI) or (hemi_p not in HEMI):
        raise ValueError("Hemispheres not found!")
    if hemi != hemi_p:
        raise ValueError("Inconsistent hemispheres in input files!")

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

    # load surfaces
    vtx_white, fac = read_geometry(surf_white)
    vtx_pial, _ = read_geometry(surf_pial)

    # load index files
    ind_white = np.loadtxt(ind_white1).astype(int)
    ind_white2 = np.loadtxt(ind_white2).astype(int)
    ind_white = ind_white[ind_white2]

    ind_pial = np.loadtxt(ind_pial1).astype(int)
    ind_pial2 = np.loadtxt(ind_pial2).astype(int)
    ind_pial = ind_pial[ind_pial2]

    # match vertex numbers
    vtx_white, vtx_pial, fac, _ = match_vertex_number(
        vtx_white, vtx_pial, fac, ind_white, ind_pial
    )

    file_white = os.path.join(path_dense, f"{hemi}.{basename_white}_match")
    file_pial = os.path.join(path_dense, f"{hemi}.{basename_pial}_match")

    write_geometry(file_white, vtx_white, fac)
    write_geometry(file_pial, vtx_pial, fac)

    # extract main component
    extract_main_component(file_white, f"{file_white}_final")
    extract_main_component(file_pial, f"{file_pial}_final")

    # inflation
    inflate_surf_mesh(
        f"{file_white}_final",
        f"{file_white}_final_inflate{n_iter_inflate}",
        n_iter_inflate,
    )

    # meshlines
    vtx_white, fac = read_geometry(f"{file_white}_final")
    vtx_pial, _ = read_geometry(f"{file_pial}_final")

    vtx_line, fac_line = get_meshlines(vtx_pial, vtx_white)
    file_line = os.path.join(path_dense, f"{hemi}.mesh_lines")
    write_geometry(file_line, vtx_line, fac_line)

    # layer
    for j in range(n_layer):
        vtx_layer = vtx_white + j / (n_layer - 1) * (vtx_pial - vtx_white)
        file_layer = os.path.join(path_layer, f"{hemi}.layer_{j}")
        write_geometry(file_layer, vtx_layer, fac)


def _main(argv=None):
    """Layer generation workflow."""
    options = _get_parser().parse_args(argv)
    kwargs = vars(options)
    layer_workflow(**kwargs)


if __name__ == "__main__":
    _main()
