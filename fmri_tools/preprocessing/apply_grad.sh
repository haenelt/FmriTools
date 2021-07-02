#!/bin/bash

# This script calls the function gradient_unwarp.py for gradient non-linearity 
# undistortion. The method is taken from the repository https://github.com/
# Washington-University/gradunwarp.git. All arguments which can be passed to 
# gradient_unwarp can be seen by calling gradient_unwarp with --help as 
# argument. Since we use high resolution data, we set numpoints=128 and interp_
# order=2. It has to be noted that gradient_unwarp is written in python 2.7 and 
# I include this function in my preprocessing pipeline which is written in 
# python >= 3.6. Therefore, I only call this function via a bash scripts which 
# first calls a python2 conda environment and sets it back to python3 after 
# running the function (I know, it is a little cumbersome). I.e., two conda 
# environments have to be defined with different python versions. The repository 
# should be installed in the python2 environment.
# Inputs:
#   *arg1: name of python3 environment
#   *arg2: name of python2 environment
#   *arg3: output directory
#   *arg4: input file
#   *arg5: output file
#   *arg6: coefficient file
# Outputs:
#   *trilinear: unwarped image.
#   *fullWarp_abs: warp field.

# get python2 environment
echo -e "\nPassed arguments"
echo -e "Arg1: "$1 # name of python3 environment
echo -e "Arg2: "$2 # name of python2 environment
echo -e "Arg3: "$3 # output directory
echo -e "Arg4: "$4 # input file
echo -e "Arg5: "$5 # output file
echo -e "Arg6: "$6 # coefficient file

# load python2 environment
source ~/.bashrc
conda activate $2

# which python version
python_version="$(python -V 2>&1)"
echo -e "\nPython version: "$python_version"\n"

# change to working directory as gradient_unwarp.py outputs some files directly into pwd
cd $3

# gradient unwarping
gradient_unwarp.py $4 $5 siemens -g $6 -n --numpoints 128 --interp_order 2

# load python3 environment
source ~/.bashrc
conda activate $1

# which python version
python_version="$(python -V 2>&1)"
echo -e "\nPython version: "$python_version
