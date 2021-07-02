% Multipol analysis
%
% This script runs a retinotopy analysis for computing the BOLD amplitude
% from a phase-encoding paradigm. The analysis includes baseline correction
% in temporal domain followed by calculating the amplitude at stimulus
% frequency in Fourier domain. If a baseline correction was already done
% (i.e. if a file with prefix b exists), no baseline correction is
% performed. If the number of cycles <freq> is not an integer, the first
% cycle fraction is discarded from analysis. Cleanup will delete all
% generated time series at the end.

% input data
input.data = {
    '/data/pt_01880/Experiment4_PSF/p6/psf/SE_EPI2/multipol_2/udata.nii',...
    '/data/pt_01880/Experiment4_PSF/p6/psf/SE_EPI2/multipol_4/udata.nii',...
    '/data/pt_01880/Experiment4_PSF/p6/psf/SE_EPI2/multipol_6/udata.nii',...
    '/data/pt_01880/Experiment4_PSF/p6/psf/SE_EPI2/multipol_8/udata.nii',...
    '/data/pt_01880/Experiment4_PSF/p6/psf/SE_EPI2/multipol_10/udata.nii',...
    '/data/pt_01880/Experiment4_PSF/p6/psf/SE_EPI2/multipol_12/udata.nii',...
    '/data/pt_01880/Experiment4_PSF/p6/psf/SE_EPI2/multipol_14/udata.nii',...
    }; % anticlock
input.tr = 3; % repetition time in s
input.period = 48; % cycle period in s
input.fix = 12; % pre and post run baseline block in s
input.freq = 10.5; % number of cycles
input.cutoff = 144; % cutoff frequency 1/cutoff in Hz

% output specification
name_sess = 'SE_EPI2';
path_output = '/data/pt_01880/Experiment4_PSF/p6/psf';
cleanup = true;

%%% do not edit below %%%

for i = 1:length(input.data)

    % path and filename
    [path, file, ext] = fileparts(input.data{i});

    % run baseline correction
    if ~exist(fullfile(path,['b' file ext]), 'file')
        baseline_correction(...
            input.data{i},...
            input.tr,...
            input.cutoff);
    end

    % run frequency analysis
    retino_fourier(...
        fullfile(path, ['b' file ext]),...
        input.freq,...
        input.fix,...
        input.period,...
        input.tr);

    % run multipol phase
    multipol_phase(...
        fullfile(path,['rb' file '_real' ext]),...
        fullfile(path,['rb' file '_imag' ext]),...
        input.freq,...
        name_sess,...
        path_output);
    
    if cleanup
        delete(fullfile(path,['b' file ext]));
        delete(fullfile(path,['rb' file '_real' ext]));
        delete(fullfile(path,['rb' file '_imag' ext]));
    end
    
end
