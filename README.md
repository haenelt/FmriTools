FmriTools
===

[![Python](https://img.shields.io/badge/Python-3.6%7C3.7%7C3.8-blue)](https://github.com/haenelt/GBB)
[![License](https://img.shields.io/badge/license-GPL--3.0-orange)](https://github.com/haenelt/GBB)

Python package (including some Matlab functions) for processing and analyzing high-resolution fMRI data. Various scripts are included which I currently use for data processing. Please be aware that all functions are written for my own convenience are are under contiuous development.

## Installation
I recommend to use `Anaconda` to create a new python environment with `Python >= 3.6`. Then, clone this repository and run the following line from the directory in which the repository was cloned with the environment being activated:

```
python setup.py install
```

Some scripts need optional (non-default) packages which can be installed with the following command:

```
pip install fmri_tools[addon]
```

Another option is to install the package in development mode by creating a `conda.pth` (which includes the path of the cloned repository) in the site-packages folder of your conda environment. Necessary external python packages can be installed by running `pip install -r requirements.txt`. With `pip install -r requirements_extra.txt`, optional (non-default) packages are installed which are only needed for single scripts.

## Other dependencies
The package cannot be import without the installation of two further modules:

`pip install pycortex`<br>
`pip install nighres`

Additionally, the installation of FreeSurfer [[1]](#1), FSL [[2]](#1), AFNI [[3]](#1), ANTS [[4]](#1), LAYNII [[5]](#1) and SPM12 [[6]](#6) is needed in some functions. 

Furthermore, these scripts call additional toolboxes:

1. `./scripts/processing/do_hmri.m` needs an installation of the hMRI toolbox [[7]](#7).
2. `./scripts/processing/vasa_calibration.m` needs an installation of the VASA toolbox [[8]](#8).
3. `./scripts/processing/TDM.m` needs installation of knkutils [[9]](#9), GLMdenoise [[10]](#10) and TDM [[11]](#11).
4. `./scripts/registration/recursiveBBR.py` needs an installation of the OpenFmriAnalysis toolbox [[12]](#12).
5. `./scripts/processing/alff.py` needs an installation of the PhysIO Toolbox [[13]](#13).
6. `./scripts/preprocessing/gnl_undistortion.py` and `./scripts/preprocessing/gnl_undistortion_time.py` need an installation of gradunwarp [[14]](#14).

## References
<a id="1">[1]</a> https://surfer.nmr.mgh.harvard.edu/<br>
<a id="2">[2]</a> https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/<br>
<a id="3">[3]</a> https://afni.nimh.nih.gov/<br>
<a id="4">[4]</a> https://www.nitrc.org/projects/ants<br>
<a id="5">[5]</a> https://github.com/layerfMRI/LAYNII<br>
<a id="6">[6]</a> https://www.fil.ion.ucl.ac.uk/spm/software/spm12/<br>
<a id="7">[7]</a> Kazan S et al. Physiological basis of vascular autocalibration (VasA): Comparison to hypercapnia calibration methods, MRM 78(3), 1168--1173 (2017).<br>
<a id="8">[8]</a> https://hmri-group.github.io/hMRI-toolbox/<br>
<a id="9">[9]</a> https://github.com/kendrickkay/knkutils<br>
<a id="11">[10]</a> https://github.com/kendrickkay/GLMdenoise<br>
<a id="12">[11]</a> https://github.com/kendrickkay/TDM<br>
<a id="13">[12]</a> https://github.com/TimVanMourik/OpenFmriAnalysis<br>
<a id="14">[13]</a> https://www.nitrc.org/projects/physio/<br>
<a id="15">[14]</a> https://github.com/Washington-University/gradunwarp<br>

## Contact
If you have questions, problems or suggestions regarding the FmriTools package, please feel free to contact [me](mailto:daniel.haenelt@gmail.com).