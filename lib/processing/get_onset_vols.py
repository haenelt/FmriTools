def get_onset_vols(cond_input, outlier_input, name_condition, TR, skip_vol):
    """
    This function returns all the volume indices corresponding to an experimental condition in from
    a block design.
    Inputs:
        *cond_input: block design condition .mat file.
        *outlier_input: regressor of no interest .txt file.
        *TR: repetition time in s.
        *skip_vol: number of skipped time point of each block.
        *name_condition: name of first condition.
    Outputs:
        *onsets: sorted volumes of experimental condition.
    
    created by Daniel Haenelt
    Date created: 16-09-2019
    Last modified: 17-03-2020
    """
    import sys
    import numpy as np
    from scipy.io import loadmat

    # load condition file
    cond = loadmat(cond_input)

    # get condition information
    names = np.concatenate(np.concatenate(cond["names"]))
    onsets = np.concatenate(cond["onsets"])
    durations = np.concatenate(np.concatenate(np.concatenate(cond["durations"])))

    # check if condition names exist
    if not name_condition in names:
        sys.exit("The condition is not found in the condition_file")

    # index of condition
    ind = np.where(name_condition == names)[0][0]

    # onsets for both conditions
    onsets = np.round(onsets[ind][0] / TR + skip_vol).astype(int)

    # durations for both conditions
    durations = np.round(durations[ind] / TR - skip_vol).astype(int)

    # sort all volumes to be considered in both conditions
    temp = onsets.copy()
    for i in range(durations - 1):
        onsets = np.append(onsets, temp + i + 1)
    onsets = np.sort(onsets)
    
    # remove outlier volumes
    if outlier_input:
        
        # load outlier regressor
        outlier_regressor = np.loadtxt(outlier_input)
        outlier_regressor = np.where(outlier_regressor == 1)[0]
        
        # look for outliers in onset arrays
        for i in range(len(onsets)):
            if np.any(onsets[i] == outlier_regressor):
                onsets[i] = -1
        
        # remove outliers
        onsets = onsets[onsets != -1]
    
    return onsets