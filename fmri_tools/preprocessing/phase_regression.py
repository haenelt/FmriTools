"""
# General
- [1]
    - post-processing identification and removal of venous signals using a phase 
      regressor technique
    - temporal magnitude and phase signal changes from larger veins in response to task 
      have an approximately linear relationship
    - such phase changes were identified and filtered from the magnitude time course, 
      resulting in suppression of macroscopic venous effects
    - the phase regressor estimates the signal in the magnitude that is explainable by 
      the phase data, S_est, by finding linear fit parameters A and B, S_est = A*phi+B
    - Fitereing is performed by subtracting out the estimate from the original data: 
      S_filt = S-S_est
    - Magnitude and phase time series in MRI tend to have different SNR characteristics 
      due to the underlying noise disitributions
    - from a data-fitting viewpoint, this results in different measurement uncertainties 
      in magnitude and phase
    - it is therfore impoartant to condition the fit of the phase regressor with 
      knowledge of the relative error levels in the measurements
    - As a coarse measure, one can estimate the standard deviations sigma_s and 
      sigma_phi of the magnitude and phase time series, respectively, if the signal 
      changes of interest can be factored out
    - use temporal high-frequency part of time series for sigma estimation
    - data fitting was performed in python using the scipy.odr interface to the Fortran
      ODRPACK library
- [2]
    - phase of the GE-EPI is sensitive to vessel size and this provides a potential 
      avenue to reduce the macrovascular weighting of the signal (phase regression)
    - IV effects
    - EV effects should be canceled
    - However, EV frequency shifts could produce a phase change in a sufficiently small 
      voxel when the symmetry assumption is violated (Vu and Gallant, 2015)
    - ODR requires inputs to estimate error ellipses prior to fitting.
    - errors estimates for magnitude and phase by high-pass filtering at 0.15 Hz (above 
      task frequency)
    - get temporal standard deviation
- [3]
    - Because of the symmetry of the EV magnetic field perturbations (cos2phi term), the 
      EV soubs ub a vixek tebd ti gave a near-zero net phase, while the intravascular 
      spins have a coherent phase shift relative to the rotating frame
    - changes in phase are likely to be observable only if vessels giving rise to S_IV 
      are not randomly oriented (as a population) within the voxel, and have a 
      sufficient blood volume for the phase change to be detectable above the noise
    - the sigma_S and sigma_phase estimates were derived from the SD of the time series 
      by notch filtering S(i) and phase(i) at the paradigm frequency and its first four 
      harmonics, and assuming the remaining variation was due to noise
    - sigma_phi and sigma_m are determined individually for each voxel by taking the FT 
      of the time series, setting the task frequency and its first four harmonics to 
      zero, and then calculating the variance of the filtered time series after taking 
      the IFT
    - it is clear that increase or decrease of BOLD signal magnitude can give rise to 
      phase changes of either sign, presumably depending on the voxel's vascular 
      orientation with respect to the field
    - it has been shown empirically in this study that regardless of the size of vessels 
      in the voxel and their orientation, a linear fit suffices
- [4]
    - source-localized phase regression (sPR) builds upon the phase regressor (PR) large 
      vein suppression method developed by Menon, 2002 with two key improvements. First, 
      sPR utilizes and unbiased least squares loss function to avoid the overcorrection 
      of large vein contributions as reported by Nencka and Rowe (2007). Second, sPR 
      utilizes task-related phase changes in neighboring voxels to suppress large vein 
      contributions even in voxels with poor phase SNR
    - However, Nencka and Rowe (2007) cautioned against the use of these large vein 
      suppression techniques given their tendency to over- or undercorrect for large 
      vein contributions. Furthermore, these large vein suppression techniques assume 
      that voxels containing large veins will have high phase fSNR. This assumption 
      begins to breakfown for large veins roughly the size of a voxel as well as for 
      veins near the magic angle
    - sPR uses unbiased ordinary least squares loss function
    - To address the poor suppression of large veins roughly the size of a voxel and 
      those near the magic angle, sPR takes advantage of the magnitude and phase fSNR 
      distbution around such veins. Voxels that contain a large vein roughly the size of 
      a voxel have high magnitude but low phase fSNR. This is because the off resonance 
      field distortions sampled by these ovxels vary symmetrically around 0 Hz. 
      Furthermore, IV task-related phase changes in these voxels, especially at higher 
      field strengths, will contirbute minimally due to the short T2* of venous blood 
      relative to that of the parenchyma.
    - sPR uses neighboring voxels (voxel i and 6 face neighbors) whose phase component 
      is most correlated wit hthe magnitude component of voxel i
    - they use two runs (one run for correlation estimation) to avoid circularity
- [5]
    - Savitzky-Golay Filtering (SGF) requires two parameters to be specified
    - N: frame size
    - p: polynomial order
    - The optimal values for N and p are unknown and vary from one voxel to the next, so 
      we implemented an exploratory implementarion of SGF that considered many possible 
      combinations of p and N for every voxel
    - We define "optimal" as the combination of N and p that minimizes the temporal 
      variance of the magnitude time series and thus produces the highest r^2 (where the 
      goodness o ffit metric r^2=1-sigma_pr/sigma_orig) where sigma_orig denotes 
      temporal standard deviation of a voxel before processing and sigm_pr denotes the 
      temporal standard deviation of the same voxel after PR
    - For EPI, we considered 11 values of p (2, 3, 4, ..., 11, 12) and up to 12 values 
      of N (5, 9, 13, ..., 45, 49) with the condition N > p
    - We do not make a priori assumpations regarding the evolution of phase in any voxel 
      because it can be modulated by both Delta R_2^* and spatiotemporally varying B0 
      inhomogeneities. We therefore adopt an aggressive data-driven approach to identify 
      the SGF parameters that minimize temporal variance (i.e., noise) for each voxel
    - By summary, we have identified that increased physiological noise at ultra-high 
      fields can significantly impede the efficacy of vessel suppression via PR, but 
      that this efficacy can be restored after low-pass filtering the phase time series 
      with a SGF
    
# Literature
[1] Curtis et al. (2014): Phase based venous suppression in resting-state BOLD GE-fMRI
[2] Stanley et al. (2021): Effects of phase regression on high-resolution functional MRI 
    of the primary visual cortex
[3] Menon (2002): Postacquisition suppression of large-vessel BOLD signals in high-
    resolution fMRI
[4] Vu and Gallant (2015): Using a novel source-localized phase regressor technique for 
    evaluation of the vascular contribution to semantic category area localization in 
    BOLD fMRI
[5] Barry and Gore (2014): Enhanced phase regression with Savitzky-Golay filtering for 
    high-resolution BOLD MRI

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
    - apply motion parameters to real and imaginary parts
    - since it is spatially smooth and interpolatable
    - convert back to phase
    - unwrap and detrend data
    - subtracting off the value of the first volume, resulting in delta-phase time-
      course
    - mask data using magnitude mask
    - calculate snr (mean, noise, snr)

3. Regression
    - PhaseFitODR.py
    - regress magnitude and phase
    - map of correlation of magnitude and phase
    - map of r^2
    - map of residuals
    - resting-state and task-based
"""
