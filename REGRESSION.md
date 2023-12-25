# Literature
- Curtis (2014): Phase based venous suppression in resting-state BOLD GE-fMRI
- Menon (2002): Postacquisition suppression of large-vessel BOLD signals in high-resolution fMRI
- Stanley (2021): Effects of phase regression on high-resolution functional MRI of the primary visual cortex

# Reference repository
- https://github.com/brettliem/phaseprep

# Procedure
1. Preprocessing of magnitude image
    - convert image to float
    - motion correct each run
    - calculate relative motions
    - generate motion plots
    - linear detrending
    - brain extraction
    - calculate snr (mean, noise, snr)

2. Preprocessing phase image
    - convert data to float
    - determine scaling required for radians
    - apply radian scaling
    - apply magnitude motion correction parameters
    - unwrap and detrend data
    - mask data using magnitude mask
    - calculate snr (mean, noise, snr)

3. Regression
    - regress magnitude and phase
    - get correlation of magnitude and phase
    - get residuals
    - resting-state and task-based

