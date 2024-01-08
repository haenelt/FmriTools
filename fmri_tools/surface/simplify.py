# -*- coding: utf-8 -*-
"""Mesh simplification
ADD
- edge collapse
- remove edges which are longer than a defined threshold
- both vertices are merged into one location at the mid-point of the removed edge
"""

import os
import subprocess
import sys

from ..io.filename import get_filename

__all__ = ["extract_main_component"]


def extract_main_component(file_in, file_out):
    """This function removes unconnected parts found in a surface mesh and returns
    the main component.

    Parameters
    ----------
    file_in : str
        Filename of input surface.
    file_out : str
        Filename of output surface.

    Returns
    -------
    None.

    """
    # make output folder
    path_output, _, _ = get_filename(file_out)
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # extract main component
    try:
        subprocess.run(["mris_extract_main_component", file_in, file_out], check=True)
    except subprocess.CalledProcessError:
        sys.exit("Main component extraction failed!")
