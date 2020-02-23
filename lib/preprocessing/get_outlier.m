function outlier_all = get_outlier(file_rp, file_in, outlier_params, path_output)
% This function uses the realignment parameter file to identify motion
% outlier volumes which exceed a defined motion threshold. Furthermore,
% the realigned time series is checked for intensity outliers by checking
% if the average z-score of each time point exceeds a defined threshold. A
% summary file with the listed outliers and a regressor of no interest are
% written.
% Inputs:
    % file_rp: filename of textfile with realignment parameters.
    % file_in: filename of realigned time series.
    % outlier_params: structure containing outlier thresholds.
    % path_output: path where output is written.
    
% created by Daniel Haenelt
% Date created: 23-02-2020
% Last modified: 23-02-2020

% make output folder
if ~exist(path_output,'dir') 
    mkdir(path_output); 
end

% initialize outlier arrays
motion_outlier = struct();
motion_outlier.short.mm = [];
motion_outlier.long.mm = [];
motion_outlier.short.deg = [];
motion_outlier.long.deg = [];
intensity_outlier = [];
    
% read realignment parameters
M = dlmread(file_rp);

% rad2deg
M(:,4) = radtodeg(M(:,4));
M(:,5) = radtodeg(M(:,5));
M(:,6) = radtodeg(M(:,6));

% translational displacement
M_trans = sqrt(M(:,1).^2+M(:,2).^2+M(:,3).^2);

% rotational displacement
M_rot = sqrt(M(:,4).^2+M(:,5).^2+M(:,6).^2);

% check for translational outliers
for i = 1:length(M_trans)-1
    diff_short = abs(M_trans(i+1) - M_trans(i));
    diff_long = abs(M_trans(i+1) - M_trans(1));

    if diff_short >= outlier_params.moco_out_mm_short
        motion_outlier.short.mm = [motion_outlier.short.mm ; [i+1 diff_short]];
    end
    
    if diff_long >= outlier_params.moco_out_mm_long
        motion_outlier.long.mm = [motion_outlier.long.mm ; [i+1 diff_long]];
    end
end

% check for rotational outliers
for i = 1:length(M_rot)-1
    diff_short = abs(M_rot(i+1) - M_rot(i));
    diff_long = abs(M_rot(i+1) - M_rot(1));

    if diff_short >= outlier_params.moco_out_deg_short
        motion_outlier.short.deg = [motion_outlier.short.deg ; [i+1 diff_short]];
    end
    
    if diff_long >= outlier_params.moco_out_deg_long
        motion_outlier.long.deg = [motion_outlier.long.deg ; [i+1 diff_long]];
    end
end

% read time series
data_img = spm_vol(file_in);
data_array = spm_read_vols(data_img);
nt = length(data_img);

data_mean = mean(data_array(:));
for i = 1:nt
    % compute threshold based on z-score
    data_array(data_array < data_mean) = NaN;
    z = (data_array(:,:,:,i) - nanmean(data_array,4)) ./ nanstd(data_array,0,4);
    z_avg = nanmean(z(:));

    if z_avg > outlier_params.int_out_z
        intensity_outlier = [intensity_outlier ; [i z_avg]];
    end
end

% get regressor of not interest
outlier_all = [];

if motion_outlier.short.mm
    outlier_all = [outlier_all motion_outlier.short.mm(:,1)];
end

if motion_outlier.long.mm
    outlier_all = [outlier_all motion_outlier.long.mm(:,1)];
end
   
if motion_outlier.short.deg
    outlier_all = [outlier_all motion_outlier.short.deg(:,1)];
end

if motion_outlier.long.deg
    outlier_all = [outlier_all motion_outlier.long.deg(:,1)];
end

if intensity_outlier
    outlier_all = [outlier_all intensity_outlier(:,1)];
end

outlier_all = unique(outlier_all);

regressor = zeros(nt,1);
regressor(outlier_all) = 1;

% get basename of output
[~, name_output, ~] = fileparts(file_in);

% write regressor of no interest
dlmwrite(fullfile(path_output,['outlier_regressor_' name_output '.txt']), regressor);

% save summed outliers in mat file
save(fullfile(path_output,['outlier_' name_output '.mat']),'motion_outlier','intensity_outlier');