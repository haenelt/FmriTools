% Preprocessing of fmri time series

% This script performs slice timing correction, fieldmap undistortion and
% realignment of an array of fmri time series. Realignment can be done
% across multiple time series. Note that the same slice timing and fieldmap 
% parameters are assumed for all time series. So, it makes most sense to 
% only put time series from the same session into one array.
%
% BandwidthPerPixelPhaseEncode: 
% retinotopy_1p0 -> 20.27
% resting_state_0p8 and task_0p8 -> 16.304
% localiser_2p0 -> 20.027

% created by Daniel Haenelt
% Date created: 06-08-2019
% Last modified: 18-03-2020

% array of of input time series
img_input = {
    '/home/daniel/Schreibtisch/test_sess/Run_1/data.nii',...
    '/home/daniel/Schreibtisch/test_sess/Run_2/data.nii',...
    };

% slice timing parameters
slice_params.slice_timing = false; % run slice timing correction
slice_params.TR = 3; % repetition time in seconds
slice_params.slice_order = 'descending'; % slice ordering (ascending or descending)

% fieldmap parameters
field_params.fieldmap_undistortion = false; % run fieldmap undistortion
field_params.fmap_magn = '';
field_params.fmap_phase = '';
field_params.fmap_te1 = 6.0; % shorter echo time in ms
field_params.fmap_te2 = 7.02; % longer echo time in ms
field_params.fmap_blipdir = -1; % phase-encoding direction
field_params.fmap_BandwidthPerPixelPhaseEncode = 16.304; % phase-encoding bandwidth Hz/px

% realignment parameters
realign_params.unwarp = false;
realign_params.mask = false;
realign_params.c = [95 130 25];
realign_params.r = [35 25 20];

% outlier parameters
outlier_params.moco_out_mm_short = 0.4; % in mm
outlier_params.moco_out_mm_long = 0.8; % in mm
outlier_params.moco_out_deg_short = 0.5; % in rad
outlier_params.moco_out_deg_long = 1.0; % in rad
outlier_params.int_out_z = 2; % in z-score

% data range parameters
range_params.apply = true;
range_params.data_min = 0;
range_params.data_max = 4095;

% add spm and lib to path
pathSPM = '/data/pt_01880/source/spm12';
pathLIB = '/data/hu_haenelt/projects/scripts/lib/preprocessing';

%%% do not edit below %%%

% add paths to the interpreter's search path
addpath(pathLIB);

% start preprocessing
fmri_preprocessing(...
    img_input,...
    slice_params,...
    field_params,...
    realign_params,...
    outlier_params,...
    range_params,...
    pathSPM);
