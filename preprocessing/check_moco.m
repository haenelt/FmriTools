% Check motion and intensity outliers

% This script recomputes the outlier estimation from fmri time series
% preprocessing. The input image array is a list of preprocessed time
% series (with the prefix u).

% created by Daniel Haenelt
% Date created: 15-09-2019
% Last modified: 27-09-2019

% array of of input time series
img_input = {
    '/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI2/Run_1/uadata.nii',....
    '/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI2/Run_2/uadata.nii',....
    '/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI2/Run_3/uadata.nii',....
    '/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI2/Run_4/uadata.nii',....
    '/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI2/Run_5/uadata.nii',....
    '/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI2/Run_6/uadata.nii',....
    '/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI2/Run_7/uadata.nii',....
    '/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI2/Run_8/uadata.nii',....
    '/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI2/Run_9/uadata.nii',....
    '/data/pt_01880/Experiment2_Rivalry/p3/rivalry/GE_EPI2/Run_10/uadata.nii',....
    };

% outlier parameters
outlier_params.moco_out_mm_short = 0.4; % in mm
outlier_params.moco_out_mm_long = 0.8; % in mm
outlier_params.moco_out_rad_short = 0.02; % in rad
outlier_params.moco_out_rad_long = 0.04; % in rad
outlier_params.int_out_z = 2; % in z-score

% add spm to path
pathSPM = '/data/pt_01880/source/spm12';

%%% do not edit below %%%

% add spm to path
addpath(pathSPM);
spm('defaults','FMRI');
spm_get_defaults('stats.maxmem',2^35);% maxmen indicates how much memory can be used
spm_get_defaults('cmdline',true);% no gui

% path diagnosis
if length(img_input) > 1
    path_diagnosis = fullfile(fileparts(fileparts(img_input{1})),'diagnosis');
else
    path_diagnosis = fullfile(fileparts(img_input{1}),'diagnosis');
end

if ~exist(path_diagnosis,'dir') 
    mkdir(path_diagnosis); 
end

% delete old preprocessing summary
[~, file, ~] = fileparts(img_input{1});
file = file(2:end); % get rid of prefix u

if exist(fullfile(path_diagnosis,['preprocessing_summary_' file '.txt']), 'file') == 2
    delete(fullfile(path_diagnosis,['preprocessing_summary_' file '.txt']));
end

% get motion and intensity outliers
for i = 1:length(img_input)

    % paths
    path_regressor = fullfile(fileparts(img_input{i}),'logfiles');
    if ~exist(path_regressor,'dir') 
        mkdir(path_regressor); 
    end
    
    % delete old outlier summary
    if exist(fullfile(fileparts(img_input{i}),['outlier_' file '.mat']), 'file') == 2
        delete(fullfile(fileparts(img_input{i}),['outlier_' file '.mat']));
    end
    
    % delete old motion regressor
    if exist(fullfile(path_regressor,['outlier_regressor_' file '.txt']), 'file') == 2
        delete(fullfile(path_regressor,['outlier_regressor_' file '.txt']));
    end
    
    % initialise outlier array
    motion_outlier = struct();
    motion_outlier.short.mm = [];
    motion_outlier.short.dir = [];
    motion_outlier.short.t = [];
    motion_outlier.long.mm = [];
    motion_outlier.long.dir = [];
    motion_outlier.long.t = [];
    
    % read realignment parameters
    M = dlmread(fullfile(fileparts(img_input{i}),['rp_' file '.txt']));
        
    % compute absolute difference between volumes
    for j = 1:length(M(:,1))-1
        for k = 1:length(M(1,:)) % number of motion parameters
            diff_short = abs(M(j+1,k) - M(j,k));
            diff_long = abs(M(j+1,k) - M(1,k));
            
            if k < 4 % set different threshold for displacement and rotation
                threshold_short = outlier_params.moco_out_mm_short;
                threshold_long = outlier_params.moco_out_mm_long;
            else
                threshold_short = outlier_params.moco_out_rad_short;
                threshold_long = outlier_params.moco_out_rad_long;
            end
            
            % check for volume-to-volume outliers
            if diff_short >= threshold_short
                motion_outlier.short.mm = [motion_outlier.short.mm diff_short];
                motion_outlier.short.dir = [motion_outlier.short.dir k];
                motion_outlier.short.t = [motion_outlier.short.t j+1];
            end
            
            % check for volume-to-reference outliers
            if diff_long >= threshold_long
                motion_outlier.long.mm = [motion_outlier.long.mm diff_long];
                motion_outlier.long.dir = [motion_outlier.long.dir k];
                motion_outlier.long.t = [motion_outlier.long.t j+1];
            end
        end
    end
    
    % initialise outlier array
    intensity_outlier = struct();
    intensity_outlier.z = [];
    intensity_outlier.t = [];
    
    % read time series
    data_img = spm_vol(img_input{i});
    data_array = spm_read_vols(data_img);
    nt = length(data_img);
    
    data_mean = mean(data_array(:));
    for j = 1:nt
        % compute threshold based on z-score
        data_array(data_array < data_mean) = NaN;
        z = (data_array(:,:,:,j) - nanmean(data_array,4)) ./ nanstd(data_array,0,4);
        z_avg = nanmean(z(:));
        
        if z_avg > outlier_params.int_out_z
            intensity_outlier.z = [intensity_outlier.z z_avg];
            intensity_outlier.t = [intensity_outlier.t j];
        end
    end
    
    % save summed outliers in mat file
    save(fullfile(fileparts(img_input{i}),['outlier_' file '.mat']),'motion_outlier','intensity_outlier');
end

% open textfile
fileID = fopen(fullfile(path_diagnosis,['preprocessing_summary_' file '.txt']),'w');
fprintf(fileID,'Percentage of within-run motion and intensity outliers\n');
fprintf(fileID,['motion threshold (mm, short): ' num2str(outlier_params.moco_out_mm_short) '\n']);
fprintf(fileID,['motion threshold (mm, long): ' num2str(outlier_params.moco_out_mm_long) '\n']);
fprintf(fileID,['motion threshold (rad, short): ' num2str(outlier_params.moco_out_rad_short) '\n']);
fprintf(fileID,['motion threshold (rad, long): ' num2str(outlier_params.moco_out_rad_long) '\n']);
fprintf(fileID,['intensity threshold (z-score): ' num2str(outlier_params.int_out_z) '\n']);
fprintf(fileID,'----------\n\n');

for i  = 1:length(img_input)

    load(fullfile(fileparts(img_input{i}),['outlier_' file '.mat']));
    nt = length(spm_vol(img_input{i}));
    path_regressor = fullfile(fileparts(img_input{i}),'logfiles');
    
    % get within-run outlier percentage
    outlier_all = [motion_outlier.short.t motion_outlier.long.t intensity_outlier.t];
    outlier_all = unique(outlier_all);
    outlier_percentage = length(outlier_all) / nt * 100;
    
    % get ratio of outliers within time series
    fprintf(fileID,'%.2f\n',outlier_percentage);

    M = zeros(nt,1);
    M(outlier_all) = 1;
    dlmwrite(fullfile(path_regressor,['outlier_regressor_' file '.txt']),M);
    
end

% close file
fclose(fileID);