% VASA estimation

% This script computes baseline corrected contrast files using the VasA
% fMRI approach for task-based data.

% created by Daniel Haenelt
% Date created: 12-02-2019
% Last modified: 12-02-2019

% input data
img_input = {
    '/data/pt_01880/MonocularBlindness/p1/flicker/GE_EPI2/Run_1/udata.nii',...
    '/data/pt_01880/MonocularBlindness/p1/flicker/GE_EPI2/Run_2/udata.nii',...
    '/data/pt_01880/MonocularBlindness/p1/flicker/GE_EPI2/Run_3/udata.nii',...
    '/data/pt_01880/MonocularBlindness/p1/flicker/GE_EPI2/Run_4/udata.nii',...
    '/data/pt_01880/MonocularBlindness/p1/flicker/GE_EPI2/Run_5/udata.nii',...
    '/data/pt_01880/MonocularBlindness/p1/flicker/GE_EPI2/Run_6/udata.nii',...
    '/data/pt_01880/MonocularBlindness/p1/flicker/GE_EPI2/Run_7/udata.nii',...
    };

contrast_input = {
    '/data/pt_01880/MonocularBlindness/p1/flicker/results/spmT/native/spmT_off_on_GE_EPI2.nii',...
    '/data/pt_01880/MonocularBlindness/p1/flicker/results/spmT/native/spmT_on_off_GE_EPI2.nii',...
    };

spm_mat_input = '/data/pt_01880/MonocularBlindness/p1/flicker/GE_EPI2/contrast/SPM.mat';

%%% do not edit below %%%

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
