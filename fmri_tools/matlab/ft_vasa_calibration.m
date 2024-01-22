function ft_vasa_calibration(img_input, contrast_input, spm_mat_input)
% VASA estimation
%
% ft_vasa_calibration(img_input, contrast_input, spm_mat_input)
% 
% Inputs:
%   img_input      - cell array of filenames of input time series.
%   contrast_input - cell array of SPM contrasts.
%   spm_mat_input  - filename of SPM.mat file.
% 
% This function computes baseline corrected contrast files using the VasA
% fMRI approach for task-based data.
%
% Dependencies (should be installed as SPM toolbox):
% - VASA toolbox: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5573956/

% add paths to the interpreter's search path
spm('defaults','FMRI');
spm_get_defaults('stats.maxmem',2^35); % maxmen indicates how much memory can be used
spm_get_defaults('cmdline',true); % no gui

% number of volumes is taken from the first entry of the input list
data_img = spm_vol(img_input{1});
nt = length(data_img);

% VasA estimation
matlabbatch{1}.spm.tools.vasa.estimate.spmmat = {spm_mat_input};
for i = 1:length(img_input)
    for j = 1:nt
        matlabbatch{1}.spm.tools.vasa.estimate.data.scans{1,i}{j,1} = [img_input{i} ',' num2str(j)];
    end
end
matlabbatch{1}.spm.tools.vasa.estimate.freq = [0.01 0.08];
matlabbatch{1}.spm.tools.vasa.estimate.fwhm = [0 0 0];

% run
spm_jobman('run',matlabbatch);

% clear matlabbatch
clear matlabbatch

% VasA apply
matlabbatch{1}.spm.tools.vasa.apply.vasamap = {fullfile(fileparts(spm_mat_input),'sRALFF_tfMRI_4D_ResidualMaps.nii,1')};
for i = 1:length(contrast_input)
    matlabbatch{1}.spm.tools.vasa.apply.contrasts{i,1} = [contrast_input{i} ',1'];
end

% run
spm_jobman('run',matlabbatch);

% clear matlabbatch
clear matlabbatch
