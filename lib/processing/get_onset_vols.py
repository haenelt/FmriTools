def get_onset_vols(cond_input, outlier_input, condition1, condition2, TR, skipvol):
    """
    This function blas...
    Inputs:
        *cond_input: block design condition .mat file.
        *outlier_input: regressor of no interest .txt file.
        *condition1: name of first condition.
        *condition2: name of second condition.
        *TR: time serier repetition time in s.
        *skipvol: number of skipped time point of each block.
    Outputs:
        *onsets1: sorted volumes of first condition.
        *onsets2: sorted volumes of second condition.
    
    created by Daniel Haenelt
    Date created: 16-09-2019
    Last modified: 16-09-2019
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
    if not condition1 in names:
        sys.exit("condition1 not found in the condition_file")

    if not condition2 in names:
        sys.exit("condition2 not found in the condition_file")

    # index of both conditions    
    c1 = np.where(condition1 == names)[0][0]
    c2 = np.where(condition2 == names)[0][0]

    # onsets for both conditions
    onsets1 = np.round(onsets[c1][0] / TR + skipvol).astype(int)
    onsets2 = np.round(onsets[c2][0] / TR + skipvol).astype(int)

    # durations for both conditions
    durations1 = np.round(durations[c1] / TR - skipvol).astype(int)
    durations2 = np.round(durations[c2] / TR - skipvol).astype(int)

    # sort all volumes to be considered in both conditions
    temp = onsets1.copy()
    for j in range(durations1 - 1):
        onsets1 = np.append(onsets1, temp + j + 1)
    onsets1 = np.sort(onsets1)

    temp = onsets2.copy()
    for j in range(durations2 - 1):
        onsets2 = np.append(onsets2, temp + j + 1)
    onsets2 = np.sort(onsets2)
    
    # remove outlier volumes
    if outlier_input:
        
        # load outlier regressor
        outlier_regressor = np.loadtxt(outlier_input)
        outlier_regressor = np.where(outlier_regressor == 1)[0]
        
        # look for outliers in onset arrays
        for j in range(len(onsets1)):
            if np.any(onsets1[j] == outlier_regressor):
                onsets1[j] = -1
        
        for j in range(len(onsets2)):
            if np.any(onsets2[j] == outlier_regressor):
                onsets2[j] = -1
        
        # remove outliers
        onsets1 = onsets1[onsets1 != -1]
        onsets2 = onsets2[onsets2 != -1]
    
    return onsets1, onsets2