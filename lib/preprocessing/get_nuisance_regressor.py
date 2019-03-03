def get_nuisance_regressor(input, wm_mask, csf_mask, path_output):
    """
    This function creates nuisance regressors from a functional time series using wm and csf masks.
    Inputs:
        *inupt: (baseline corrected) time series.
        *wm_mask: white matter mask registered to the time series.
        *csf_mask: csf mask registered to the time series.
        *path_output: path where output is saved.

    created by Daniel Haenelt
    Date created: 02-03-2019    
    Last modified: 02-03-2019
    """
    import os
    import numpy as np
    import nibabel as nb
    import matplotlib.pyplot as plt

    # make output folder
    if not os.path.exists(path_output):
        os.mkdir(path_output)
    
    # get mean signal in mask
    wm_array = nb.load(wm_mask).get_fdata()
    csf_array = nb.load(csf_mask).get_fdata()

    # make regressor
    func_array = nb.load(input).get_fdata()
    nt = np.shape(func_array)[3]

    # get ROI mean signal
    nuisance_regressor = np.zeros((nt,2))
    for i in range(nt):
        nuisance_regressor[i,0] = np.mean(func_array[:,:,:,i][wm_array == 1])
        nuisance_regressor[i,1] = np.mean(func_array[:,:,:,i][csf_array == 1])

    # save regressor
    np.savetxt(os.path.join(path_output,"nuisance_regressor.txt"),nuisance_regressor,'%.7e','\t')

    # plot regressor
    plt.figure(1, figsize=(12, 6))
    plt.plot(nuisance_regressor[:,0])
    plt.title('White matter variation')
    plt.ylabel('Intensity in a.u.')
    plt.xlabel('Volume')
    plt.xticks(np.arange(0, len(nuisance_regressor[:,0])))
    plt.savefig(os.path.join(path_output,'wm_regressor.png'))

    plt.figure(2, figsize=(12, 6))
    plt.plot(nuisance_regressor[:,1])
    plt.title('CSF variation')
    plt.ylabel('Intensity in a.u.')
    plt.xlabel('Volume')
    plt.xticks(np.arange(0, len(nuisance_regressor[:,1])))
    plt.savefig(os.path.join(path_output,'csf_regressor.png'))