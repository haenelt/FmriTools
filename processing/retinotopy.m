% Retinotopy analysis

% This script runs a retinotopy analysis for a phase encoding paradigm
% with rings and wedges. The analysis includes the following steps: 
% (1) baseline correction, (2) phase calculation, (3) averaging over 
% sessions with opposite stimulus direction. If a baseline correction was
% already done (i.e. if a file with prefix b exists), no baseline
% correction is performed. Note that the number of cycles <freq> is not an
% integer value in my runs. The first cycle fraction is discarded from
% further analysis. We use this to have already initialised activation in
% the first volume.

% created by Daniel Haenelt
% Date created: 07-12-2018
% Last modified: 18-02-2019

% input data
input.pol.data.pos = '/data/pt_01880/Experiment3_Stripes/p3/retinotopy2/pol_anticlock/udata_gnlcorr.nii'; % anticlock
input.pol.data.neg = '/data/pt_01880/Experiment3_Stripes/p3/retinotopy2/pol_clock/udata_gnlcorr.nii'; % clock
input.pol.tr = 2; % repetition time in s
input.pol.period = 64; % cycle period in s
input.pol.fix = 12; % pre run baseline block in s (post run baseline did not work)
input.pol.freq = 8.25; % number of cycles
input.pol.cutoff = 192; % cutoff frequency 1/cutoff in Hz

input.ecc.data.pos = '/data/pt_01880/Experiment3_Stripes/p3/retinotopy2/ecc_expanding/udata_gnlcorr.nii'; % expanding
input.ecc.data.neg = '/data/pt_01880/Experiment3_Stripes/p3/retinotopy2/ecc_contracting/udata_gnlcorr.nii'; % contracting
input.ecc.tr = 2;
input.ecc.period = 32;
input.ecc.fix = 12;
input.ecc.freq = 8.25;
input.ecc.cutoff = 96;

% add spm to path
pathSPM = '/data/pt_01880/source/spm12'; 

% add library to path
pathLIB = '/home/raid2/haenelt/projects/scripts/lib/';

%%% do not edit below %%%

% add paths to the interpreter's search path
addpath(pathSPM);
cd(fileparts(which(mfilename)));
addpath(fullfile(pathLIB,'preprocessing'));
addpath(fullfile(pathLIB,'processing'));

% concatenate all runs into single cell
all_data = {input.pol.data.pos input.pol.data.neg input.ecc.data.pos input.ecc.data.neg};
all_tr = [input.pol.tr input.pol.tr input.ecc.tr input.ecc.tr];
all_cutoff = [input.pol.cutoff input.pol.cutoff input.ecc.cutoff input.ecc.cutoff];
all_freq = [input.pol.freq input.pol.freq input.ecc.freq input.ecc.freq];
all_fix = [input.pol.fix input.pol.fix input.ecc.fix input.ecc.fix];
all_period = [input.pol.period input.pol.period input.ecc.period input.ecc.period];

% prepare path and filename
path = cell(1,length(all_data));
file = cell(1,length(all_data));
ext = cell(1,length(all_data));
for i = 1:length(all_data)
    [path{i}, file{i}, ext{i}] = fileparts(all_data{i});
end

% run baseline correction
for i = 1:length(all_data)    
    if exist(fullfile(path{i},['b' file{i} ext{i}]), 'file')
        continue
    else
        baseline_correction(...
            all_data{i},...
            all_tr(i),...
            all_cutoff(i),...
            pathSPM);
    end
end

% run frequency analysis
for i = 1:length(all_data)
    retino_fourier(...
        fullfile(path{i}, ['b' file{i} ext{i}]),...
        all_freq(i),...
        all_fix(i),...
        all_period(i),...
        all_tr(i),...
        pathSPM);
end

% run phase averaging
[pos_path, pos_file, pos_ext] = fileparts(input.pol.data.pos);
[neg_path, neg_file, neg_ext] = fileparts(input.pol.data.neg);
average_phase(...
    fullfile(pos_path,['rb' pos_file '_real' pos_ext]),...
    fullfile(pos_path,['rb' pos_file '_imag' pos_ext]),...
    fullfile(neg_path,['rb' neg_file '_real' neg_ext]),...
    fullfile(neg_path,['rb' neg_file '_imag' neg_ext]),...  
    'pol',...
    input.pol.freq,...
    pathSPM);

[pos_path, pos_file, pos_ext] = fileparts(input.ecc.data.pos);
[neg_path, neg_file, neg_ext] = fileparts(input.ecc.data.neg);
average_phase(...
    fullfile(pos_path,['rb' pos_file '_real' pos_ext]),...
    fullfile(pos_path,['rb' pos_file '_imag' pos_ext]),...
    fullfile(neg_path,['rb' neg_file '_real' neg_ext]),...
    fullfile(neg_path,['rb' neg_file '_imag' neg_ext]),...  
    'ecc',...
    input.ecc.freq,...
    pathSPM);