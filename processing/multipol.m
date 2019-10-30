% Multipol analysis

% This script runs a retinotopy analysis for computing the BOLD amplitude
% from a phase-encoding paradigm. The analysis includes baseline correction
% in temporal domain followed by calculating the amplitude at stimulus
% frequency in Fourier domain. If a baseline correction was already done
% (i.e. if a file with prefix b exists), no baseline correction is
% performed. If the number of cycles <freq> is not an integer, the first
% cycle fraction is discarded from analysis.

% created by Daniel Haenelt
% Date created: 04-03-2019
% Last modified: 12-07-2019

% input data
input.data = {
    '/data/pt_01880/Experiment4_PSF/p6/psf/GE_EPI1/multipol_2/udata.nii',...
    '/data/pt_01880/Experiment4_PSF/p6/psf/GE_EPI1/multipol_4/udata.nii',...
    '/data/pt_01880/Experiment4_PSF/p6/psf/GE_EPI1/multipol_6/udata.nii',...
    '/data/pt_01880/Experiment4_PSF/p6/psf/GE_EPI1/multipol_8/udata.nii',...
    '/data/pt_01880/Experiment4_PSF/p6/psf/GE_EPI1/multipol_10/udata.nii',...
    '/data/pt_01880/Experiment4_PSF/p6/psf/GE_EPI1/multipol_12/udata.nii',...
    '/data/pt_01880/Experiment4_PSF/p6/psf/GE_EPI1/multipol_14/udata.nii',...
    }; % anticlock
input.tr = 3; % repetition time in s
input.period = 48; % cycle period in s
input.fix = 12; % pre run baseline block in s
input.freq = 10.5; % number of cycles
input.cutoff = 144; % cutoff frequency 1/cutoff in Hz

% output specification
name_sess = 'GE_EPI1';
path_output = '/data/pt_01880/Experiment4_PSF/p6/psf';

% add spm and lib to path
pathSPM = '/data/pt_01880/source/spm12'; 
pathLIB = '/home/raid2/haenelt/projects/scripts/lib';

%%% do not edit below %%%

% add paths to the interpreter's search path
addpath(pathSPM);
addpath(fullfile(pathLIB,'preprocessing'));
addpath(fullfile(pathLIB,'processing'));

for i = 1:length(input.data)

    % path and filename
    [path, file, ext] = fileparts(input.data{i});

    % run baseline correction
    if ~exist(fullfile(path,['b' file ext]), 'file')
        baseline_correction(...
            input.data{i},...
            input.tr,...
            input.cutoff,...
            pathSPM);
    end

    % run frequency analysis
    retino_fourier(...
        fullfile(path, ['b' file ext]),...
        input.freq,...
        input.fix,...
        input.period,...
        input.tr,...
        pathSPM);

    % run multipol phase
    multipol_phase(...
        fullfile(path,['rb' file '_real' ext]),...
        fullfile(path,['rb' file '_imag' ext]),...
        input.freq,...
        pathSPM,...
        name_sess,...
        path_output);
    
end