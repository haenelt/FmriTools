"""
Gradient nonlinearity correction

This scripts calls the HCP toolbox to correct for gradient nonlinearities in the input volume.

Before running the script, login to queen via ssh and set the fsl environment by calling FSL in 
the terminal.

created by Daniel Haenelt
Date created: 10-01-2020
Last modified: 10-01-2020
"""
from lib.preprocessing.gnl_correction import gnl_correction

# input
input = [
    "/data/pt_01880/Experiment3_Stripes/p1/mpm/pd_kp_mtflash3d_v1ax_0p5_0008/Results/s2899875-124739-00001-00352-1_PD.nii",
    "/data/pt_01880/Experiment3_Stripes/p1/mpm/pd_kp_mtflash3d_v1ax_0p5_0008/Results/s2899875-124739-00001-00352-1_R1.nii",
    "/data/pt_01880/Experiment3_Stripes/p1/mpm/pd_kp_mtflash3d_v1ax_0p5_0008/Results/s2899875-124739-00001-00352-1_R2s_WOLS.nii",
    ]

path_output = [
    "/data/pt_01880/Experiment3_Stripes/p1/mpm/pd_kp_mtflash3d_v1ax_0p5_0008/Results",
    ]

file_bash = "/home/raid2/haenelt/projects/gradunwarp/apply_grad.sh"
file_coeff = "/home/raid2/haenelt/projects/gradunwarp/7t_coeff.grad"
python3_env = "daniel"
python2_env = "daniel2"
cleanup = True

""" do not edit below """

for i in range(len(input)):
    gnl_correction(input[i], file_bash, file_coeff, python3_env, python2_env, path_output[i], 
                   cleanup)
