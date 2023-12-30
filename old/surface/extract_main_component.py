# -*- coding: utf-8 -*-

import os
import subprocess
import sys

from ..io.filename import get_filename


def extract_main_component(file_in, file_out):
    """Extract main component.

    This function removes unconnected parts found in a surface mesh and returns
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
