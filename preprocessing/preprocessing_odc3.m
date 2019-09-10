% Preprocessing of fmri time series

% This script performs slice timing correction, fieldmap undistortion and
% realignment of an array of fmri time series. Realignment can be done
% across multiple time series or for each time series independently. Note
% that in both cases, the same slice timing and fieldmap parameters are
% assumed for all time series. So, it makes most sense to only put time
% series from the same session into one array.
%
% BandwidthPerPixelPhaseEncode: 
% retinotopy_1p0 -> 20.27
% resting_state_0p8 and task_0p8 -> 16.304
% localiser_2p0 -> 20.027

% created by Daniel Haenelt
% Date created: 06-08-2019
% Last modified: 06-09-2019

% array of of input time series
img_input = {
    '/data/pt_01880/Experiment2_Rivalry/p3/odc/GE_EPI3/Run_1/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/odc/GE_EPI3/Run_2/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/odc/GE_EPI3/Run_3/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/odc/GE_EPI3/Run_4/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/odc/GE_EPI3/Run_5/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/odc/GE_EPI3/Run_6/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/odc/GE_EPI3/Run_7/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/odc/GE_EPI3/Run_8/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/odc/GE_EPI3/Run_9/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/odc/GE_EPI3/Run_10/data.nii',...
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
field_params.fmap_blipdir = 1; % phase-encoding direction
field_params.fmap_BandwidthPerPixelPhaseEncode = 16.304; % phase-encoding bandwidth Hz/px

% realignment parameters
realign_params.independent = false;

% outlier parameters
outlier_params.moco_out_mm_short = 0.8; % in mm
outlier_params.moco_out_mm_long = 1.6; % in mm
outlier_params.moco_out_rad_short = 0.1; % in rad
outlier_params.moco_out_rad_long = 0.2; % in rad
outlier_params.int_out_z = 2; % in z-score

% add spm and lib to path
pathSPM = '/data/pt_01880/source/spm12';
pathLIB = '/home/raid2/haenelt/projects/scripts/lib/preprocessing';

%%% do not edit below %%%

% add paths to the interpreter's search path
addpath(pathLIB);

% start preprocessing
if realign_params.independent 
    for i = 1:length(img_input)
        fmri_preprocessing(...
            {img_input{i}},...
            slice_params,...
            field_params,...
            outlier_params,...
            pathSPM);
    end
else
    fmri_preprocessing(...
        img_input,...
        slice_params,...
        field_params,...
        outlier_params,...
        pathSPM);
end
