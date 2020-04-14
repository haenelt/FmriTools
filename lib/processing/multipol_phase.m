function multipol_phase(input_real, input_imag, freq, pathSPM, name_sess, path_output)

% This function takes the real and imaginary parts from the time series
% fft. Similar to Dumoulin et al., 2017, we compute an coherence estimate
% of the voxel respose at stimulus frequency by taking the amplitude and 
% dividing it by the root mean square of the one-sided power spectrum.
% Other measures as F-statistics are computed by taking the amplitude and 
% dividing it by the mean amplitude of the rest of the frequency spectrum
% Based on this, the voxel population can be thresholded to investigate
% different subpopulations. To run this function, the number of runs has to
% be specified to choose an appropriate results folder location.
% Inputs:
    % input_real: real part of fft.
    % input_imag: imaginary part of.
    % freq: number of cycles.
    % pathSPM: path to spm toolbox.
    % name_sess: name of sessiincluded in the basename.
    % path_output: path where output is saved.

% created by Daniel Haenelt
% Date created: 04-03-2019
% Last modified: 12-07-2019

% add spm to path
addpath(pathSPM);

% parameters
% freq_ignored is a vector containing the frequencies that should be
% ignored when calculating the coherence values.
freq = floor(freq);
freq_ignored = [];
%freq_ignored = [0:1 freq-1 freq freq+1];

% define output path
[~, name_run, ~] = fileparts(fileparts(input_real));
path_output = fullfile(path_output,'results','native');
if ~exist(path_output,'dir') 
    mkdir(path_output); 
end

% load input
data_img = spm_vol(input_imag); % imaginary
data_imag_array = spm_read_vols(data_img);
data_img = spm_vol(input_real); % real
data_real_array = spm_read_vols(data_img);

% get image dimension
nt = length(data_img);

% mean power over all (considered) frequencies
pw = sqrt(data_real_array.^2+data_imag_array.^2); 
mfs = 1:round(nt/2);
mfs = mfs(~ismember(mfs, freq_ignored+1));
sumpw = sqrt(nansum(pw(:,:,:,mfs).^2,4));
mpw = nanmean(pw(:,:,:,mfs),4); 
mstd = nanstd(pw(:,:,:,mfs),0,4);

% real and imaginary parts at stimulus frequency
data_imag_freq_array = squeeze(data_imag_array(:,:,:,freq+1));
data_real_freq_array = squeeze(data_real_array(:,:,:,freq+1));   

% data at cycle frequency
A_array = sqrt(data_real_freq_array.^2+data_imag_freq_array.^2); % amplitude
F_array = A_array ./ mpw * 100; % F-statistic 
snr_array = A_array ./ mstd; % SNR
coherence_array = A_array ./ sumpw;

% transform to phases
pha_array = atan2(data_imag_freq_array,data_real_freq_array)/pi*180; % phase in degrees
pha_array = mod(pha_array,360); % phases between 0 and 360
pha_array = pha_array - 180;  % phases between -180 and +180

% change NaN in background to 0
data_real_freq_array(isnan(data_real_freq_array)) = 0;
data_imag_freq_array(isnan(data_imag_freq_array)) = 0;
pha_array(isnan(pha_array)) = 0;
A_array(isnan(A_array)) = 0;
F_array(isnan(F_array)) = 0;
snr_array(isnan(snr_array)) = 0;
coherence_array(isnan(coherence_array)) = 0;

% save niftis
nhdr = data_img(1);

if isempty(name_sess)
    nhdr.fname = fullfile(path_output,['real_' name_run '.nii']);
    spm_write_vol(nhdr,data_real_freq_array);

    nhdr.fname = fullfile(path_output,['imag_' name_run '.nii']);
    spm_write_vol(nhdr,data_imag_freq_array);

    nhdr.fname = fullfile(path_output,['phase_' name_run '.nii']);
    spm_write_vol(nhdr,pha_array);

    nhdr.fname = fullfile(path_output,['a_' name_run '.nii']);
    spm_write_vol(nhdr,A_array);

    nhdr.fname = fullfile(path_output,['f_' name_run '.nii']);
    spm_write_vol(nhdr,F_array);

    nhdr.fname = fullfile(path_output,['snr_' name_run '.nii']);
    spm_write_vol(nhdr,snr_array);

    nhdr.fname = fullfile(path_output,['c_' name_run '.nii']);
    spm_write_vol(nhdr,coherence_array);
else
    nhdr.fname = fullfile(path_output,['real_' name_sess '_' name_run '.nii']);
    spm_write_vol(nhdr,data_real_freq_array);

    nhdr.fname = fullfile(path_output,['imag_' name_sess '_' name_run '.nii']);
    spm_write_vol(nhdr,data_imag_freq_array);

    nhdr.fname = fullfile(path_output,['phase_' name_sess '_' name_run '.nii']);
    spm_write_vol(nhdr,pha_array);

    nhdr.fname = fullfile(path_output,['a_' name_sess '_' name_run '.nii']);
    spm_write_vol(nhdr,A_array);

    nhdr.fname = fullfile(path_output,['f_' name_sess '_' name_run '.nii']);
    spm_write_vol(nhdr,F_array);

    nhdr.fname = fullfile(path_output,['snr_' name_sess '_' name_run '.nii']);
    spm_write_vol(nhdr,snr_array);

    nhdr.fname = fullfile(path_output,['c_' name_sess '_' name_run '.nii']);
    spm_write_vol(nhdr,coherence_array);
end