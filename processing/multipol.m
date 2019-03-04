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
% Last modified: 04-03-2019

% input data
input.data = '/nobackup/actinium1/haenelt/V2STRIPES/p6/retinotopy/pol_anticlock/udata.nii'; % anticlock
input.tr = 3; % repetition time in s
input.period = 60; % cycle period in s
input.fix = 12; % pre run baseline block in s (post run baseline did not work)
input.freq = 8.25; % number of cycles
input.cutoff = 120; % cutoff frequency 1/cutoff in Hz

% add spm to path
pathSPM = '/nobackup/actinium1/haenelt/source/spm12'; 

% add library to path
pathLIB = '/home/raid2/haenelt/scripts/lib';

%%% do not edit below %%%

% add paths to the interpreter's search path
addpath(pathSPM);
addpath(fullfile(pathLIB,'preprocessing'));
addpath(fullfile(pathLIB,'processing'));

% path and filename
[path, file, ext] = fileparts(input.data);

% run baseline correction
if ~exist(fullfile(path,['b' file ext]), 'file')
    baseline_correction(...
        input.data,...
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
    pathSPM);