def plot_moco(path_moco):
    """
    Plot motion parameters and maximum displacement (AFNI parameters). The time course of motion
    parameters is saved as image file in the same directory.
    Inputs:
        *path_moco: path to AFNI motion parameter files.
    
    created by Daniel Haenelt
    Date created: 21-11-2018             
    Last modified: 21-11-2018
    """
    import os
    import numpy as np
    import matplotlib.pyplot as plt

    # load motion parameters
    moco_params = np.loadtxt(os.path.join(path_moco,'moco_params.1D'))

    roll = moco_params[:,0]
    pitch = moco_params[:,1]
    yaw = moco_params[:,2]
    superior = moco_params[:,3]
    left = moco_params[:,4]
    posterior = moco_params[:,5]

    # load maximum displacement parameters
    max_disp = np.loadtxt(os.path.join(path_moco,'max_disp.1D'))

    # load displacement change parameters
    delta_disp = np.loadtxt(os.path.join(path_moco,'max_disp_delt.1D'))

    # plot motion parameters
    plt.figure(1, figsize=(12, 6))
    plt.plot(roll)
    plt.plot(pitch)
    plt.plot(yaw)
    plt.title('Head rotation')
    plt.ylabel('Displacement in Degrees')
    plt.xlabel('Volume')
    plt.xticks(np.arange(0, len(roll)))
    plt.legend(['$\Delta$ roll','$\Delta$ pitch','$\Delta$ yaw'])
    plt.savefig(os.path.join(path_moco,'head_rotation.png'))
    
    plt.figure(2, figsize=(12, 6))
    plt.plot(superior)
    plt.plot(left)
    plt.plot(posterior)
    plt.title('Head translation')
    plt.ylabel('Displacement in mm')
    plt.xlabel('Volume')
    plt.xticks(np.arange(0, len(superior)))
    plt.legend(['$\Delta$ Superior','$\Delta$ Left','$\Delta$ Posterior'])
    plt.savefig(os.path.join(path_moco,'head_translation.png'))
    
    plt.figure(3, figsize=(12, 6))
    plt.plot(max_disp)
    plt.title('Maximal displacement')
    plt.ylabel('Displacement in mm')
    plt.xlabel('Volume')
    plt.xticks(np.arange(0, len(max_disp)))
    plt.savefig(os.path.join(path_moco,'max_disp.png'))
    
    plt.figure(4, figsize=(12, 6))
    plt.plot(delta_disp)
    plt.title('Displacement change between neighbouring volumes')
    plt.ylabel('$\Delta$ Displacement in mm')
    plt.xlabel('Volume')
    plt.xticks(np.arange(0, len(delta_disp)))
    plt.savefig(os.path.join(path_moco,'delta_disp.png'))
    
    plt.close('all')