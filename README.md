FmriTools
===

[![Python](https://img.shields.io/badge/Python-3.6%2B-blue)](https://github.com/haenelt/FmriTools)
[![License](https://img.shields.io/github/license/haenelt/FmriTools)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/172737126.svg)](https://zenodo.org/doi/10.5281/zenodo.10573041)

Python package for processing and analyzing high-resolution fMRI data. Please be aware that all functions are written for my own convenience and are under continuous development.

## Installation
I recommend to use `Miniconda` to create a new python environment with `Python >= 3.6`. Then, clone this repository and run the following line from the directory in which the repository was cloned with the environment being activated:

```
pip install .
```

Another option is to install the package in development mode by creating a `conda.pth` (which includes the path of the cloned repository) in the site-packages folder of your conda environment. Necessary external python packages can then be installed by running `pip install -r requirements.txt`.

The package contains some matlab functions. Please add the following code to the matlab `startup.m` script. This includes the package to the matlab search path. Next to the root directory of this package, you should add the root directory of [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/), which is used throughout the package. Of course you have to specify the paths there.

```matlab
%------------ FmriTools ------------------------------%

% Root directory of the fmri_tools package
path_fmri = <PATH-TO-FMRITOOLS-FOLDER>;

% Root directory of the SPM12 toolbox
path_spm12 = <PATH-TO-SPM12-ROOT>;

% add paths
addpath(fullfile(path_fmri, 'matlab'));
addpath(genpath(path_spm12));

clear path_fmri path_spm12;
%-----------------------------------------------------%
```

Single functions also call the matlab toolboxes [OpenFmriAnalysis](https://github.com/TimVanMourik/OpenFmriAnalysis), [knkutils](https://github.com/cvnlab/knkutils), [GLMdenoise](https://github.com/cvnlab/GLMdenoise), [TDM](https://github.com/cvnlab/TDM) and the SPM12 toolboxes [hMRI](https://hmri-group.github.io/hMRI-toolbox/), [PhysIO](https://www.tnu.ethz.ch/en/software/tapas/documentations/physio-toolbox), and [VASA](https://pubmed.ncbi.nlm.nih.gov/26416648/).

Furthermore, some functions need an installation of [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/), [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/), [AFNI](https://afni.nimh.nih.gov/), [ANTS](https://www.nitrc.org/projects/ants) and/or [LAYNII](https://github.com/layerfMRI/LAYNII) which need to be accessible from the command line.

## How I installed gradunwarp
In the following, I list all necessary steps to install the [gradunwarp toolbox](https://github.com/Washington-University/gradunwarp) to run gradient nonlinearity corrections with the scripts `gnl_undistortion.py` or `gnl_undistortion_time.py`:

1. clone the gradunwarp repository `git clone https://github.com/Washington-University/gradunwarp.git`.
2. create a second conda environment with `python=2.7`: `conda create -n <env> python=2.7`.
3. activate that environment: `conda activate <env>`.
4. install necessary packages: `pip install -r requirements.txt`. The requirements file is found in the same folder as this readme.
5. change to the cloned gradunwarp repository.
6. install the repository with `python setup.py install` while the conda environment is still activated.
7. `wget "https://github.com/Washington-University/HCPpipelines/raw/stable/global/scripts/GradientDistortionUnwarp.sh"`.
- N.B. call the module `gradient_unwarp.py` with arguments `--numpoints 128` and `--interp_order 2` to have good accuracy when working with high resolution data.

## Contact
If you have questions, problems or suggestions regarding the FmriTools package, please feel free to contact [me](mailto:daniel.haenelt@gmail.com).
