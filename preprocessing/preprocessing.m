% Preprocessing of fmri time series
%
% This script performs slice time correction, fieldmap undistortion and
% motion correction in the SPM12 framework which can be applied to a
% session consisting of multiple runs. Slice time correction and fieldmap 
% undistortion are optional.
% Fieldmap undistortion is performed without skullstripping. The total
% readout time is taken as the inverse of the BandwidthPerPixelPhaseEncode
% in Hz/px which can be found in the dicom tag (0019, 1028).
% From the realignment procedure, outliers are detected based on motion and
% intensity variations. A volume is classified if either volume-to-volume
% motion (short) or volume-to-reference motion (long) exceeds a defined
% threshold. The intensity outlier classification is taken from Chaimow et
% al., 2016. For every volume within a time-series, the fMRI signals in
% single voxels are compared to the entire time-course of that voxel by
% computing the z-score of the measured fMRI signal relative to the
% voxels's time series. If the average z-score across all voxels within one
% volume exceeds a defined threshold, this volume is marked as outlier. 
% In the outlier summary, the number of outliers within one run is compared
% to its length is computed and if the ratio exceeds a defined threshold, 
% the whole run is discarded. Additionally, regressors for neglecting 
% outlier volumes are created.
% 
% BandwidthPerPixelPhaseEncode: 
% retinotopy_1p0 -> 20.27
% resting_state_0p8 and task_0p8 -> 16.304
% localiser_2p0 -> 20.027
%
% created by Daniel Haenelt
% Date created: 26-02-2019
% Last modified: 04-07-2019

% input time series
img_input = {
    '/data/pt_01880/Experiment2_Rivalry/p4/odc/GE_EPI1/Run_1/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p4/odc/GE_EPI1/Run_2/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p4/odc/GE_EPI1/Run_3/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p4/odc/GE_EPI1/Run_4/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p4/odc/GE_EPI1/Run_5/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p4/odc/GE_EPI1/Run_6/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p4/odc/GE_EPI1/Run_7/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p4/odc/GE_EPI1/Run_8/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p4/odc/GE_EPI1/Run_9/data.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p4/odc/GE_EPI1/Run_10/data.nii',...
    };

% input fieldmap
fmap_magn = '';
fmap_phase = '';

% path to spm toolbox
pathSPM = '/data/pt_01880/source/spm12';

% slice timing parameters
slice_timing = false; % run slice timing correction
TR = 3; % repetition time in s
slice_order = 'descending'; % slice ordering (ascending or descending)

% fieldmap parameters
fieldmap_undistortion = false; % run fieldmap undistortion
fmap_te1 = 6.0; % shorter echo time in ms
fmap_te2 = 7.02; % longer echo time in ms
fmap_blipdir = 1; % phase-encoding direction
fmap_BandwidthPerPixelPhaseEncode = 16.304; % phase-encoding bandwidth Hz/px

% outlier parameters
moco_out_mm_short = 0.8; % in mm
moco_out_mm_long = 1.6; % in mm
moco_out_rad_short = 0.1; % in rad
moco_out_rad_long = 0.2; % in rad
int_out_z = 2; % in z-score
run_out = 0.1; % in percent of time series length

%%% do not edit below %%%

% start spm
addpath(pathSPM);
spm('defaults','FMRI');
spm_get_defaults('stats.maxmem',2^35);% maxmen indicates how much memory can be used
spm_get_defaults('cmdline',true);% no gui

% preprocessing summary
if length(img_input) > 1
    path_diagnosis = fullfile(fileparts(fileparts(img_input{1})),'diagnosis');
else
    path_diagnosis = fullfile(fileparts(img_input{1}),'diagnosis');
end

if ~exist(path_diagnosis,'dir') 
    mkdir(path_diagnosis); 
end

% path and filename of fieldmap phase
[path_fmap2, file_fmap2, ~] = fileparts(fmap_phase);

if slice_timing
    for i = 1:length(img_input)
        
        % length of time series
        data_img = spm_vol(img_input{i});
        nt = length(data_img); % number of volumes
        nslices = data_img(1).dim(3); % number of slices
        
        for j = 1:nt
            matlabbatch2{1}.spm.temporal.st.scans{j,1} = [img_input{i} ',' num2str(j)];
        end
        matlabbatch{1}.spm.temporal.st.nslices = nslices;
        matlabbatch{1}.spm.temporal.st.tr = TR;
        matlabbatch{1}.spm.temporal.st.ta = TR - TR / nslices;
        if strcmp(slice_order,'ascending')
            matlabbatch{1}.spm.temporal.st.so = 1:1:nslices;
            matlabbatch{1}.spm.temporal.st.refslice = 1;
        elseif strcmp(slice_order,'descending')
            matlabbatch{1}.spm.temporal.st.so = nslices:-1:1;
            matlabbatch{1}.spm.temporal.st.refslice = nslices;
        else
            disp('Only ascending and descending slice ordering implemented!');
        end
        matlabbatch{1}.spm.temporal.st.prefix = 'a';
        
        % run
        spm_jobman('run',matlabbatch);
    
        % clear
        clear matlabbatch
        
    end
end

if fieldmap_undistortion
       
    % total readout time
    total_readout = 1000/fmap_BandwidthPerPixelPhaseEncode; % in ms
     
    % calculate fieldmap
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = {fmap_phase};
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = {fmap_magn};
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [fmap_te1 fmap_te2];
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = fmap_blipdir;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = total_readout;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
    for i = 1:length(img_input)
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(i).epi = {[img_input{i} ',1']};
    end
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
    if length(img_input) > 1
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
    else
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = '';
    end
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
    
    % run
    spm_jobman('run',matlabbatch);
    
    % clear
    clear matlabbatch

end

% realignment
for i = 1:length(img_input)
    
    % get time series length
    data_img = spm_vol(img_input{i});
    nt = length(data_img);
    
    % get time series path and filename
    [path, file, ext] = fileparts(img_input{i});
    if slice_timing
        file = ['a' file];
    end
    
    % get input data
    for j = 1:nt
        matlabbatch{1}.spm.spatial.realignunwarp.data(i).scans{j,1} = fullfile(path,[file ext ',' num2str(j)]);
        if fieldmap_undistortion
            if length(img_input) > 1
                matlabbatch{1}.spm.spatial.realignunwarp.data(i).pmscan = {fullfile(path_fmap2,['vdm5_sc' file_fmap2 '_session' num2str(i) '.nii,1'])};
            else
                matlabbatch{1}.spm.spatial.realignunwarp.data(i).pmscan = {fullfile(path_fmap2,['vdm5_sc' file_fmap2 '.nii,1'])};
            end
        else
            matlabbatch{1}.spm.spatial.realignunwarp.data(i).pmscan = '';
        end
    end
    
end

% change to first run directory
cd(fileparts(img_input{1}));

% realignment (coregister only)
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 1;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 1;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 1;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 7;
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';

% estimate unwarping parameters
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';

% write unwarped images
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 7;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

% run
spm_jobman('run',matlabbatch);

% clear
clear matlabbatch

% plot motion regressors
for i = 1:length(img_input)

    % change to single run
    cd(fileparts(img_input{i}));
    
    [path, file, ext] = fileparts(img_input{i});
    if slice_timing
        file = ['a' file];
    end
    
    % read realignment parameters
    M = dlmread(['rp_' file '.txt']);
    
    transFig = figure('visible','off');
    plot(M(:,1));
    hold on
    plot(M(:,2));
    plot(M(:,3));
    title(['Translational movement in session ' num2str(i)]);
    xlabel('number of volume');
    ylabel('Translation in mm');
    legend('x','y','z');
    saveas(gcf,fullfile(path_diagnosis,['moco_mm_' file '_' num2str(i) '.png']));
    close(transFig);
    
    radFig = figure('visible','off');
    plot(M(:,4));
    hold on
    plot(M(:,5));
    hold on
    plot(M(:,6));
    hold on
    title(['Rotational movement in session ' num2str(i)]);
    xlabel('number of volume');
    ylabel('Rotation in rad');
    legend('pitch','roll','yaw');
    saveas(gcf,fullfile(path_diagnosis,['moco_rad_' file '_' num2str(i) '.png']));
    close(radFig);

end

% get motion and intensity outliers
for i = 1:length(img_input)

    % initialise outlier array
    motion_outlier = struct();
    motion_outlier.short.mm = [];
    motion_outlier.short.dir = [];
    motion_outlier.short.t = [];
    motion_outlier.long.mm = [];
    motion_outlier.long.dir = [];
    motion_outlier.long.t = [];
    
    % change to single run
    cd(fileparts(img_input{i}));
    
    [path, file, ext] = fileparts(img_input{i});
    if slice_timing
        file = ['a' file];
    end
    
    % read realignment parameters
    M = dlmread(['rp_' file '.txt']);
        
    % compute absolute difference between volumes
    for j = 1:length(M(:,1))-1
        for k = 1:length(M(1,:)) % number of motion parameters
            diff_short = abs(M(j+1,k) - M(j,k));
            diff_long = abs(M(j+1,k) - M(1,k));
            
            if k < 4 % set different threshold for displacement and rotation
                threshold_short = moco_out_mm_short;
                threshold_long = moco_out_mm_long;
            else
                threshold_short = moco_out_rad_short;
                threshold_long = moco_out_rad_long;
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
    data_img = spm_vol(fullfile(path,['u' file ext]));
    data_array = spm_read_vols(data_img);
    nt = length(data_img);
    
    data_mean = mean(data_array(:));
    for j = 1:nt
        % compute threshold based on z-score
        data_array(data_array < data_mean) = NaN;
        z = (data_array(:,:,:,j) - nanmean(data_array,4)) ./ nanstd(data_array,0,4);
        z_avg = nanmean(z(:));
        
        if z_avg > int_out_z
            intensity_outlier.z = [intensity_outlier.z z_avg];
            intensity_outlier.t = [intensity_outlier.t j];
        end
    end
    
    % save summed outliers in mat file
    save(['outlier_' file '.mat'],'motion_outlier','intensity_outlier');
end

% open textfile
fileID = fopen(fullfile(path_diagnosis,['preprocessing_summary_' file '.txt']),'w');

for i  = 1:length(img_input)

    cd(fileparts(img_input{i}));
    load(['outlier_' file '.mat']);
    
    data_img = spm_vol(img_input{i});
    nt = length(data_img);
    
    % get all within-run outliers
    outlier_all = [motion_outlier.short.t motion_outlier.long.t intensity_outlier.t];
    outlier_all = unique(outlier_all);
    
    % get ratio of outliers within time series
    outlier_percentage = length(outlier_all) / nt;
    
    if outlier_percentage >= run_out
        fprintf(fileID,'%d %d\n',length(outlier_all),1);
    else
        fprintf(fileID,'%d %d\n',length(outlier_all),0);
    end
    
    [~, file, ~] = fileparts(img_input{i});
    if slice_timing
        file = ['a' file];
    end
    
    M = dlmread(['rp_' file '.txt']);
    M = [M zeros(nt,1)];
    M(outlier_all,end) = 1;
    dlmwrite(['frp_' file '.txt'],M)
    
end

% close file
fclose(fileID);

% get time series of first volumes
data_img_out = spm_vol(img_input{1});
for i = 1:length(img_input)

    [path, file, ext] = fileparts(img_input{i});
    if slice_timing
        file = ['a' file];
    end
    
    data_img = spm_vol(fullfile(path,['u' file ext]));
    data_array = spm_read_vols(data_img);
    
    data_img_out(i).fname = fullfile(path_diagnosis, ['vol1_' file '.nii']);
    spm_write_vol(data_img_out(i), data_array(:,:,:,1));

end

% get mean data
data_img_out = spm_vol(img_input{1});
data_img_out = data_img_out(1);
data_img_out.private.dat.dim(4) = [];
data_mean_array = zeros(data_img_out(1).dim);
counter = 0;
for i = 1:length(img_input)

    [path, file, ext] = fileparts(img_input{i});
    if slice_timing
        file = ['a' file];
    end
    
    data_img = spm_vol(fullfile(path,['u' file ext]));
    data_array = spm_read_vols(data_img);
    
    data_mean_array = data_mean_array + sum(data_array,4);
    counter = counter + length(data_img);
    
end
data_mean_array = data_mean_array ./ counter;
data_img_out.fname = fullfile(path_diagnosis, ['mean_' file '.nii']);
spm_write_vol(data_img_out, data_mean_array);











