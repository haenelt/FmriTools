function multipol_noise(input_real, input_imag, freq, name_sess, ...
    path_output)
% Multipol noise
%
% multipol_noise(input_real, input_imag, freq, name_sess, ...
%    path_output)
%
% Inputs:
%   input_real  - real part of fft.
%   input_imag  - imaginary part of.
%   freq        - number of cycles.
%   name_sess   - name of session included in the basename.
%   path_output - path where output is saved.
%
% This function takes the real and imaginary parts from the time series
% fft and estimates the baseline value (noise floor) for the metrics
% computed in multipol_phase.

% created by Daniel Haenelt
% Date created: 09-03-2021
% Last modified: 09-03-2021

% parameters
% freq_ignored is a vector containing the frequencies that should be
% ignored when calculating the denominator.
freq = floor(freq);
freq_ignored = [];
%freq_ignored = [0:1 freq-1 freq freq+1];

% define output paths
[~, name_run, ~] = fileparts(fileparts(input_real));
path_noise_a = fullfile(path_output,'results','noise_a','native');
path_noise_c = fullfile(path_output,'results','noise_c','native');
path_noise_f = fullfile(path_output,'results','noise_f','native');
path_noise_snr = fullfile(path_output,'results','noise_snr','native');

if ~exist(path_noise_a,'dir') 
    mkdir(path_noise_a); 
end

if ~exist(path_noise_c,'dir') 
    mkdir(path_noise_c); 
end

if ~exist(path_noise_f,'dir') 
    mkdir(path_noise_f); 
end

if ~exist(path_noise_snr,'dir') 
    mkdir(path_noise_snr); 
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

% real and imaginary parts of all frequencies discarding the stimulus
% frequency and its first harmonic
mfs = 1:round(nt/2);
mfs = mfs(~ismember(mfs, [freq-1 freq freq+1 2*freq-1 2*freq 2*freq+1]+1));

data_imag_freq_array = nanmean(data_imag_array(:,:,:,mfs),4);
data_real_freq_array = nanmean(data_real_array(:,:,:,mfs),4);

% data at cycle frequency
A_array = sqrt(data_real_freq_array.^2+data_imag_freq_array.^2); % amplitude
F_array = A_array ./ mpw * 100; % F-statistic 
snr_array = A_array ./ mstd; % SNR
coherence_array = A_array ./ sumpw; % coherence

% change NaN in background to 0
A_array(isnan(A_array)) = 0;
F_array(isnan(F_array)) = 0;
snr_array(isnan(snr_array)) = 0;
coherence_array(isnan(coherence_array)) = 0;

% save niftis
nhdr = data_img(1);

if isempty(name_sess)
    nhdr.fname = fullfile(path_noise_a,['noise_a_' name_run '.nii']);
    spm_write_vol(nhdr,A_array);
    
    nhdr.fname = fullfile(path_noise_c,['noise_c_' name_run '.nii']);
    spm_write_vol(nhdr,coherence_array);
    
    nhdr.fname = fullfile(path_noise_f,['noise_f_' name_run '.nii']);
    spm_write_vol(nhdr,F_array);
    
    nhdr.fname = fullfile(path_noise_snr,['noise_snr_' name_run '.nii']);
    spm_write_vol(nhdr,snr_array);
else
    nhdr.fname = fullfile(path_noise_a,['noise_a_' name_sess '_' name_run '.nii']);
    spm_write_vol(nhdr,A_array);
    
    nhdr.fname = fullfile(path_noise_c,['noise_c_' name_sess '_' name_run '.nii']);
    spm_write_vol(nhdr,coherence_array);
    
    nhdr.fname = fullfile(path_noise_f,['noise_f_' name_sess '_' name_run '.nii']);
    spm_write_vol(nhdr,F_array);
    
    nhdr.fname = fullfile(path_noise_snr,['noise_snr_' name_sess '_' name_run '.nii']);
    spm_write_vol(nhdr,snr_array);
end
