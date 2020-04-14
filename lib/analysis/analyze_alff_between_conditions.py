def analyze_alff_between_conditions(input_label, input_contrast1, input_contrast2, input_rest, 
                                    min_contrast, nvert):
    """
    Comparison of resting-state between different stripes populations. Mask stripes from an input 
    contrast within a label ROI and randomly select <nvert> vertices within the final mask. An 
    independent samples t-test is computed (or Welch's test if Levene's test is significant).
    Inputs:
        *input_label: input label file to define the region of interest.
        *input_contrast1: first input contrast data for masking.
        *input_contrast2: second input contrast data for masking.
        *input_rest: input resting-state data (alff or falff).
        *min_contrast: minimum contrast for masking (t-score).
        *min_rest: minimum resting-state fluctuation within the mask.
        *nvert: number of selected vertices.
    Outputs:
        *rest1: resting-state data within condition1.
        *rest2: resting-state data within condition2.
        *t: t-score from independent samples t-test.
        *p: p-value from independent samples t-test.
        *p_levene: p-value from Levene's test.

    created by Daniel Haenelt
    Date created: 11-03-2019
    Last modified: 12-03-2019
    """
    import numpy as np
    import nibabel as nb
    import random
    from nibabel.freesurfer.io import read_label
    from scipy.stats import ttest_ind, levene

    # load data
    label = read_label(input_label)
    contrast1 = np.squeeze(nb.load(input_contrast1).get_fdata())
    contrast2 = np.squeeze(nb.load(input_contrast2).get_fdata())
    rest = np.squeeze(nb.load(input_rest).get_fdata())

    # mask data
    contrast1[contrast1 < min_contrast] = np.NaN
    contrast2[contrast2 < min_contrast] = np.NaN

    contrast1[~label] = np.NaN
    contrast2[~label] = np.NaN
    contrast1[~np.isnan(contrast1)] = 1
    contrast2[~np.isnan(contrast2)] = 1

    # get resting-state data in mask
    rest1 = rest * contrast1
    rest1 = rest1[~np.isnan(rest1)]    
    rest1 = rest1[rest1 != np.min(rest1)]

    rest2 = rest * contrast2
    rest2 = rest2[~np.isnan(rest2)]
    rest2 = rest2[rest2 != np.min(rest2)]
       
    # select random number of vertices
    label_shuffled1 = random.sample(range(0,len(rest1)), nvert)
    label_shuffled2 = random.sample(range(0,len(rest2)), nvert)

    rest1 = rest1[label_shuffled1]
    rest2 = rest2[label_shuffled2]

    # independent samples t-test
    # Levene's test is run to check for equal variances. If variances are not equal, Welch's t-test
    # is performed.
    _, p_levene = levene(rest1, rest2)
    if p_levene < 0.05:
        t, p = ttest_ind(rest1, rest2, equal_var=False)
    else:
        t, p = ttest_ind(rest1, rest2, equal_var=True)

    return rest1, rest2, t, p, p_levene