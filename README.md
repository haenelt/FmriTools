FmriTools
===

[![Python](https://img.shields.io/badge/Python-3.6%7C3.7%7C3.8-blue)](https://github.com/haenelt/GBB)
[![License](https://img.shields.io/badge/license-GPL--3.0-orange)](https://github.com/haenelt/GBB)

Various Matlab and Python scripts for processing and analyzing high-resolution fMRI data. Please be aware that these functions are written for my own convenience are are under contiuous development.

## Installation
I recommend to use `Anaconda` to create a new python environment with `Python >= 3.6`. Then, clone this repository and run the following line from the directory in which the repository was cloned with the environment being activated:

```
python setup.py install
```

Some scripts need optional (non-default) packages which can be installed with the following command:

```
python setup.py fmri_tools[addon]
```

## Other dependencies

- https://github.com/Washington-University/gradunwarp.git
- OpenFmriAnalysis
- hmritoolbox
- physio
- ants
- afni
- freesurfer
- fsl
- nighres
- laynii
- spm12
- matlab
- vasa
- knkutils (kendrick kay)
- glmdenoise (kendrick kay)
- tdm (kendrick kay)

## References

## Contact
If you have questions, problems or suggestions regarding regarding the FmriTools package, please feel free to contact [me](mailto:daniel.haenelt@gmail.com).