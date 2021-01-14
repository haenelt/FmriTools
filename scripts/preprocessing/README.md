### How I installed gradunwarp
In the following, I list all necessary steps to install the [gradunwarp toolbox](https://github.com/Washington-University/gradunwarp) to run gradient nonlinearity corrections with the scripts `gnl_undistortion.py` or `gnl_undistortion_time.py`:

1. clone the gradunwarp repository `git clone https://github.com/Washington-University/gradunwarp.git`.
2. create a second conda environment with `python=2.7`: `conda create -n <env> python=2.7`.
3. activate that environment: `conda activate <env>`.
4. install necessary packages: `pip install -r requirements.txt`. The requirements file is found in the same folder as this readme.
5. change to the cloned gradunwarp repository.
6. install the repository with `python setup.py install` while the conda environment is still activated.
7. `wget "https://github.com/Washington-University/HCPpipelines/raw/stable/global/scripts/GradientDistortionUnwarp.sh"`.
- N.B. call the module `gradient_unwarp.py` with arguments `--numpoints 128` and `--interp_order 2` to have good accuracy when working with high resolution data.