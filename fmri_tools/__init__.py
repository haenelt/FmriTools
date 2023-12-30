# -*- coding: utf-8 -*-
"""Python package for the analysis of high-resolution fMRI data."""

import os

# meta infos
__author__ = "Daniel Haenelt"
__license__ = "GPL v3"
__version__ = "2.0.0-alpha"
__status__ = "Development"

__all__ = ["SoftwareName"]


class SoftwareName:
    """Names of third-party dependencies to look up in PATH variable."""

    FREESURFER = "freesurfer"
    FSL = "fsl"
    AFNI = "afni"
    ANTS = "ants"
    LAYNII = "laynii"
    MATLAB = "matlab"


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
