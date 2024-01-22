# -*- coding: utf-8 -*-
"""Python package for the analysis of high-resolution fMRI data."""

import os
import subprocess
import sys

# meta infos
__author__ = "Daniel Haenelt"
__license__ = "GPL v3"
__version__ = "2.0.0-alpha"
__status__ = "Development"

__all__ = ["SoftwareName", "check_installation"]


class SoftwareName:
    """Names of third-party dependencies to look up in PATH variable."""

    FREESURFER = "freesurfer"
    FSL = "fsl"
    AFNI = "afni"
    ANTS = "ants"
    LAYNII = "laynii"
    MATLAB = "matlab"


def check_software():
    """Check if third-party software is installed."""
    print(70 * "-")
    print("FmriTools " + "(v" + str(__version__) + ")")
    print("author: " + str(__author__))
    print(70 * "-")

    software = SoftwareName()
    members = [
        attr
        for attr in dir(software)
        if not callable(getattr(software, attr)) and not attr.startswith("__")
    ]

    print("\nCheck installations:")
    print(20 * "-")
    for member in members:
        name = SoftwareName.__dict__[member]
        if os.environ["PATH"].find(name) != -1:
            print(f"{name:<14}: found".upper())
        else:
            print(f"{name:<14}: not found".upper())


def check_installation(command):
    """Check if command can be executed from the command line."""
    try:
        subprocess.run([command], shell=True, check=False)
    except FileNotFoundError:
        sys.exit(
            f"\nCould not find '{command}'. Make sure all required software is \
                installed and can be executed from the command line."
        )
