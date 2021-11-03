FmriTools
===

[![Python](https://img.shields.io/badge/Python-3.6%7C3.7%7C3.8-blue)](https://github.com/haenelt/FmriTools)
[![License](https://img.shields.io/github/license/haenelt/FmriTools)](https://www.gnu.org/licenses/gpl-3.0)

Python package for processing and analyzing high-resolution fMRI data. Various scripts are included which I currently use for data processing. Please be aware that all functions are written for my own convenience and are under continuous development.

## Installation
I recommend to use `Miniconda` to create a new python environment with `Python >= 3.6`. Then, clone this repository and run the following line from the directory in which the repository was cloned with the environment being activated:

```
python setup.py install
```

Some scripts need optional (non-default) packages which can be installed with the following command:

```
pip install fmri_tools[addon]
```

Another option is to install the package in development mode by creating a `conda.pth` (which includes the path of the cloned repository) in the site-packages folder of your conda environment. Necessary external python packages can be installed by running `pip install -r requirements.txt`. With `pip install -r requirements_extra.txt`, optional (non-default) packages are installed which are only needed for single scripts.

The package cannot be import without the installation of two further modules:

`pip install pycortex`<br>
`pip install nighres`

The package contains some matlab functions. Please add the following code to the matlab `startup.m` script. This includes the package to the matlab search path. Of course you have to specify the paths there. Next to the root directory of this package, you should add the root directories of [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/), [knkutils](https://github.com/kendrickkay/knkutils), [GLMdenoise](https://github.com/kendrickkay/GLMdenoise), [TDM](https://github.com/kendrickkay/TDM) and [OpenFmriAnalysis](https://github.com/TimVanMourik/OpenFmriAnalysis). SPM12 is used throughout the package. The other toolboxes are used by single scripts and might be ignored.

```matlab
%------------ FmriTools ------------------------------%

% Root directory of the fmri_tools package
path_fmri = <PATH-TO-FMRITOOLS-FOLDER>;

% Root directory of the SPM12 toolbox
path_spm12 = <PATH-TO-SPM12-ROOT>;

% Root directories of some toolboxes which are used
% sporadically. Please check the readme for further
% information.
path_knkutils = <PATH-TO-KNKUTILS-FOLDER>;
path_glmdenoise = <PATH-TO-GLMDENOISE-FOLDER>;
path_tdm = <PATH-TO-TDM-FOLDER>;
path_bbr = <PATH-TO-OPENFMRIANALYSIS-FOLDER>;

% add fmri_tools
addpath(fullfile(path_fmri, 'io'));
addpath(fullfile(path_fmri, 'preprocessing'));
addpath(fullfile(path_fmri, 'processing'));
addpath(fullfile(path_fmri, 'registration'));
addpath(fullfile(path_fmri, 'segmentation'));
addpath(fullfile(path_fmri, 'skullstrip'));

% add spm12
addpath(genpath(path_spm12));

% add further toolboxes
addpath(genpath(path_knkutils));
addpath(genpath(path_glmdenoise));
addpath(genpath(path_tdm));
addpath(genpath(path_bbr));

clear path_fmri path_spm12 path_knkutils ...
    path_glmdenoise path_tdm path_bbr;
%-----------------------------------------------------%
```

Additionally, [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/), [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/), [AFNI](https://afni.nimh.nih.gov/), [ANTS](https://www.nitrc.org/projects/ants) and [LAYNII](https://github.com/layerfMRI/LAYNII) need to be installed and included to the system path.

## Dependencies of single scripts
- The script `./scripts/processing/TDM.m` needs access to the toolboxes [knkutils](https://github.com/kendrickkay/knkutils), [GLMdenoise](https://github.com/kendrickkay/GLMdenoise), [TDM](https://github.com/kendrickkay/TDM). These toolboxes should be added to the matlab search path as explained above.
- The script `./scripts/registration/recursiveBBR.py` needs access to the [OpenFmriAnalysis](https://github.com/TimVanMourik/OpenFmriAnalysis) toolbox. This toolbox should be added to the matlab search path as explained above.
- The script `./scripts/processing/do_hmri.m` needs an installation of the [hMRI toolbox](https://hmri-group.github.io/hMRI-toolbox/). This toolbox needs to be installed as SPM toolbox.
- The script `./scripts/processing/vasa_calibration.m` needs an installation of the [VASA toolbox](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5573956/). This toolbox needs to be installed as SPM toolbox.
- The script `./scripts/processing/alff.py` needs an installation of the [PhysIO toolbox](https://www.nitrc.org/projects/physio/). This toolbox needs to be installed as SPM toolbox.
6. The scripts `./scripts/preprocessing/gnl_undistortion.py` and `./scripts/preprocessing/gnl_undistortion_time.py` need an installation of [gradunwarp](https://github.com/Washington-University/gradunwarp). You can find a summary of how I installed the toolbox [here](./scripts/preprocessing/README.md).

## Contact
If you have questions, problems or suggestions regarding the FmriTools package, please feel free to contact [me](mailto:daniel.haenelt@gmail.com).
