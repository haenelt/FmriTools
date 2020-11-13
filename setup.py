# -*- coding: utf-8 -*-

# external inputs
import setuptools


INSTALL_REQUIREMENTS = ['numpy',
                        'nibabel',
                        'scipy',
                        'matplotlib',
                        'scikit-image',
                        'imageio',
                        'shapely',
                        'descartes',
                        'nipype',
                        'surfdist',
                        'sh',
                        'h5py',
                        'gbb',
                        ]

EXTRA_REQUIREMENTS = {"addon": ['pydicom', 'natsort', 'joblib']}

CLASSIFIERS = ["Programming Language :: Python :: 3",
               "Programming Language :: Python :: 3.6",
               "Programming Language :: Python :: 3.7",
               "Programming Language :: Python :: 3.8",
               "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
               "Operating System :: OS Independent",
               "Development Status :: 3 - Alpha",
               "Intended Audience :: Science/Research",
               "Topic :: Scientific/Engineering",
               ]

with open("VERSION", "r", encoding="utf8") as fh:
    VERSION = fh.read().strip()

with open("README.md", "r", encoding="utf8") as fh:
    LONG_DESCRIPTION = fh.read()

setuptools.setup(
    name="fmri_tools",
    version=VERSION,
    author="Daniel Haenelt",
    author_email="daniel.haenelt@gmail.com",
    description="Various functions for analysis of high-resolution fMRI data",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/haenelt/FmriTools",
    license='GPL v3',
    packages=setuptools.find_packages(),
    install_requires=INSTALL_REQUIREMENTS,
    extras_require=EXTRA_REQUIREMENTS,
    classifiers=CLASSIFIERS,
    python_requires='>=3.6',
    include_package_data=True,
    zip_safe=False,
    )
