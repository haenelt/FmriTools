function ft_block_contrast(img_input, cond_input, TR, name_output, name_sess, ...
output_folder, null_input, multi_input, cutoff_highpass, microtime_onset, hrf_cbv, ...
hrf_derivative, lowpass, cutoff_lowpass, order_lowpass)
% GLM analysis
%
% ft_block_contrast(img_input, cond_input, TR, name_output, name_sessoutput_folder, ...
%    null_input, multi_input, cutoff_highpass, microtime_onset, hrf_cbv, ...
%    hrf_derivative, lowpass, cutoff_lowpass, order_lowpass)
%
% Inputs:
%   img_input       - cell array of filenames of input time series.
%   cond_input      - cell array of condition files in *.mat format.
%   TR              - repetition time in s.
%   cutoff_highpass - 1/cutoff_highpass frequency in Hz.
%   name_output     - basename of output contrasts. 
%   name_sess       - name of session (if multiple sessions exist).
%   output_folder   - name of folder where SPM-mat is saved.
%   null_input      - cell array of nulling regressors.
%   multi_input     - cell array of regressors of no interest.
%   microtime_onset - only change 1 (default: 8) if reference slice in 
%                     slice timing is first slice.
%   hrf_cbv         - use an HRF designed for CBV responses.
%   hrf_derivative  - include HRF derivative in model.
%   lowpass         - lowpass filter fMRI time series.
%   cutoff_lowpass  - cutoff frequency for lowpass filtering in Hz.
%   order_lowpass   - order of lowpass filter.
%
% This function computes a fixed-effects GLM for one or multiple runs within 
% one session. The script is designed for paradigms with 2-4 experimental 
% conditions. However, for more than 2 conditions, not all contrasts are 
% computed:) Stimulus parameters are read from a spm12 compatible *.mat 
% file. Optionally, regressors of no interest and scan nulling regressors
% can be specified. Regressors of no interest are expected to be in a
% format like the rp*.txt file. Scan nulling regressors are expected to be
% loaded as textfile containing one column with either zeroes for valid 
% time points or ones for outliers in different rows.

% default parameter
if ~exist('null_input','var')  
    null_input = {};
end

if ~exist('multi_input','var')  
    multi_input = {};
end

if ~exist('microtime_onset','var')  
    microtime_onset = 8;
end

if ~exist('hrf_cbv','var')  
    hrf_cbv = false;
end

if ~exist('hrf_derivative','var')  
    hrf_derivative = false;
end

if ~exist('lowpass','var')  
    lowpass = false;
    cutoff_lowpass = 10;
    order_lowpass = 1;
end

% add paths to the interpreter's search path
spm('defaults','FMRI');
spm_get_defaults('stats.maxmem',2^35); % maxmen indicates how much memory can be used
spm_get_defaults('cmdline',true); % no gui

% lowpass filtering
if lowpass
    for i = 1:length(img_input)
        ft_lpfilter(img_input{i}, TR, cutoff_lowpass, order_lowpass);
    end
    
    % change input to lowpass filtered time series
    [filepath, name, ext] = fileparts(img_input{i});
    img_input{i} = fullfile(filepath,['l' name ext]);
end

% output folder is taken from the first entry of the input list
if length(img_input) > 1
    path_output = fullfile(fileparts(fileparts(img_input{1})),output_folder);
else
    path_output = fullfile(fileparts(img_input{1}),output_folder);
end

if ~exist(path_output,'dir')
    mkdir(path_output);
end
cd(path_output);

% initialize vector for number of noise regressors per run
n_noise = zeros(length(img_input),1);

% number of volumes is taken from the first entry of the input list
data_img = spm_vol(img_input{1});
nt = length(data_img);

% change hrf for cbv response (personal discussion with Denis Chaimow)
if hrf_cbv == true
    spm_get_defaults('stats.fmri.hrf', [5.5, 16, 1, 1, -6, 0, 32]);
end

% fmri model
matlabbatch{1}.spm.stats.fmri_spec.dir = {path_output};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = microtime_onset;

for i = 1:length(img_input)
    
    % add data
    for j = 1:nt
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans{j,1} = [img_input{i} ',' num2str(j)];
    end
    
    % add condition regressors
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi = cond_input(i);
    
    % add scan nulling regressors
    if ~isempty(null_input)
        outlier = dlmread(null_input{i});
        outlier = find(outlier == 1);
        n_noise(i) = n_noise(i) + length(outlier); 
        if isempty(outlier)
            matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress = struct('name', {}, 'val', {});
        else
            for j = 1:length(outlier)
                scan_null_regressor = zeros(nt,1);
                scan_null_regressor(outlier(j)) = 1;
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress(j).name = ['Scan nulling regressor ' num2str(j)];
                matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress(j).val = scan_null_regressor;
            end
        end
    else
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress = struct('name', {}, 'val', {});
    end
    
    % add other regressors of no interest
    if ~isempty(multi_input)
        outlier = dlmread(multi_input{i});
        n_noise(i) = n_noise(i) + size(outlier,2);
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = multi_input(i);
    else
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = {''};
    end
    
    % add highpass filter size
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf = cutoff_highpass;    
end

% hrf model
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
if hrf_derivative == true
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1 0];
else
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
end
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.3; % -Inf for no implicit mask
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

% run
spm_jobman('run',matlabbatch);

% clear matlabbatch
clear matlabbatch

% model estimation
matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(path_output,'SPM.mat')};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

% run
spm_jobman('run',matlabbatch);

% clear after completion
clear matlabbatch

% calculate contrasts
ft_tcontrast(cond_input, path_output, name_output, name_sess, hrf_derivative, n_noise);
