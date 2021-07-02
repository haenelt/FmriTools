% Preprocessing of fmri time series
%
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

% array of of input time series
img_input = {
    '/data/pt_01880/temp_p5/Run_1/data.nii',...
    '/data/pt_01880/temp_p5/Run_2/data.nii',...
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
outlier_params.moco_out_deg_short = 0.5; % in deg
outlier_params.moco_out_deg_long = 1.0; % in deg
outlier_params.int_out_z = 2; % in z-score

% separate realignment for each time series
run_separate = false;

%%% do not edit below %%%

% start preprocessing
if run_separate 
    for i = 1:length(img_input)
        fmri_preprocessing(...
            img_input(i),...
            slice_params,...
            field_params,...
            realign_params,...
            outlier_params);
    end
else
    fmri_preprocessing(...
        img_input,...
        slice_params,...
        field_params,...
        realign_params,...
        outlier_params);
end
