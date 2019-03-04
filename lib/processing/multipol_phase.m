function multipol_phase(input_real, input_imag, freq, pathSPM)

% This function takes the real and imaginary parts from the time series
% fft. Similar to Dumoulin et al., 2017, we compute an approximate
% F-statistics by taking the amplitude of the stimulation frequency and
% dividing it by the mean amplitude of the rest of the frequency spectrum
% excluding the offset (zero frequency) and the stimulation frequency.
% Based on this, the voxel population can be thresholded to investigate
% different subpopulations.
% Inputs:
    % input_real: real part of fft.
    % input_imag: imaginary part of.
    % freq: number of cycles.
    % pathSPM: path to spm toolbox.

% created by Daniel Haenelt
% Date created: 04-03-2019
% Last modified: 04-03-2019

% add spm to path
addpath(pathSPM);

% parameters
% freq_ignored is a vector containing the frequencies that should be
% ignored when calculating the coherence values.
freq = floor(freq);
freq_ignored = [0:1 freq-1 freq freq+1];

% define output path
path_output = fileparts(input_real);

% load input
data_img = spm_vol(input_imag); % imaginary
data_imag_array = spm_read_vols(data_img);
data_img = spm_vol(input_real); % real
data_real_array = spm_read_vols(data_img);

% get image dimension
nt = length(data_img);

% mean power over all (considered) frequencies
pw = sqrt(data_real_array.^2+data_imag_array.^2); 
mfs = 1:nt;
mfs = mfs(~ismember(mfs, freq_ignored+1));
mpw = nanmean(pw(:,:,:,mfs),4); 
mstd = nanstd(pw(:,:,:,mfs),0,4);

% real and imaginary parts at stimulus frequency
data_imag_freq_array = squeeze(data_imag_array(:,:,:,freq+1));
data_real_freq_array = squeeze(data_real_array(:,:,:,freq+1));   

% data at cycle frequency
A_array = sqrt(data_real_freq_array.^2+data_imag_freq_array.^2); % amplitude
F_array = A_array ./ mpw * 100; % F-statistic 
snr_array = A_array ./ mstd; % SNR

% transform to phases
pha_array = atan2(data_imag_freq_array,data_real_freq_array)/pi*180; % phase in degrees
pha_array = mod(pha_array,360); % phases between 0 and 360
pha_array = pha_array - 180;  % phases between -180 and +180

% change NaN in background to 0
data_real_freq_array(isnan(data_real_freq_array)) = 0;
data_imag_freq_array(isnan(data_imag_freq_array)) = 0;
pha_array(isnan(avg_pha)) = 0;
A_array(isnan(A_array)) = 0;
F_array(isnan(F_array)) = 0;
snr_array(isnan(snr_array)) = 0;

% save niftis
nhdr = data_img(1);
nhdr.fname = fullfile(path_output,'real.nii');
spm_write_vol(nhdr,data_real_freq_array);

nhdr.fname = fullfile(path_output,'imag.nii');
spm_write_vol(nhdr,data_imag_freq_array);

nhdr.fname = fullfile(path_output,'phase.nii');
spm_write_vol(nhdr,pha_array);

nhdr.fname = fullfile(path_output,'a.nii');
spm_write_vol(nhdr,A_array);

nhdr.fname = fullfile(path_output,'f.nii');
spm_write_vol(nhdr,F_array);

nhdr.fname = fullfile(path_output,'snr.nii');
spm_write_vol(nhdr,snr_array);