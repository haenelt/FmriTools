function regress_physio(input, regressor, TR, cutoff_highpass, path_output, cleanup)
% This function computes nuisance regression of resting-state data from
% predefined regressors. The output is the residual time series from the
% GLM.
% Inputs:
    % input: input time series
    % regressors: multiple regressors as separate columns in a textfile.
    % TR: repetition time in s.
    % cutoff_highpass: highpass filter cutoff frequency in 1/Hz
    % path_output: path where output is saved.
    % cleanup: delete intermediate files.

% created by Daniel Haenelt
% Date created: 28-02-2019
% Last modified: 28-02-2019

% default parameter
if ~exist('cleanup','var')  
    cleanup = true;
end

% add paths to the interpreter's search path
spm('defaults','FMRI');
spm_get_defaults('stats.maxmem',2^35); % maxmen indicates how much memory can be used
spm_get_defaults('cmdline',true); % no gui

% make output folder
if ~exist(path_output,'dir') 
    mkdir(path_output); 
end

% get paths and filenames of input time series
[~, file,~] = fileparts(input);

% number of volumes
data_img = spm_vol(input);
nt = length(data_img);

% fmri model
matlabbatch{1}.spm.stats.fmri_spec.dir = {path_output};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

for i = 1:nt
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans{i,1} = [input ',' num2str(i)];
end

matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {regressor};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = cutoff_highpass;

% hrf model
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.0;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

% run
spm_jobman('run',matlabbatch);

% clear matlabbatch
clear matlabbatch

% model estimation
matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(path_output,'SPM.mat')};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 1;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

% run
spm_jobman('run',matlabbatch);

% clear after completion
clear matlabbatch

% get cleaned time series
for i = 1:nt
    res_img = spm_vol(fullfile(path_output,['Res_' num2str(i,'%04d') '.nii']));
    res_array = spm_read_vols(res_img);
    res_array(isnan(res_array)) = 0;
    
    data_img(i).dt = [16 0]; % float
    data_img(i).fname = fullfile(path_output,['r' file '.nii']);
    spm_write_vol(data_img(i),res_array);
end

% clean folder
if cleanup
    cd(path_output);
    delete('beta*.nii');
    delete('Res_*.nii');
    delete('mask.nii');
    delete('ResMS.nii');
    delete('RPV.nii');
end