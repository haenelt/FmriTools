function ft_fmri_preprocessing(img_input, slice_params, field_params, ...
    realign_params, outlier_params)
% FMRI preprocessing
%
% ft_fmri_preprocessing(img_input, slice_params, field_params, ...
%    realign_params, outlier_params)
%
% Inputs:
%   input          - cell array of filenames of input time series.
%   slice_params   - struct of slice time correction parameters.
%   field_params   - struct of fieldmap parameters.
%   realign_params - struct of realignment parameters.
%   outlier_params - struct of realignment check parameters.
%
% This function performs slice time correction, fieldmap undistortion and
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
% In the outlier summary, the percentage of within-run outliers are printed
% out. Additionally, regressors of no interest for neglecting outlier
% volumes are created.
%
% Example input
%
% array of of input time series
%img_input = {
%    '/data/pt_01880/temp_p5/Run_1/data.nii',...
%    '/data/pt_01880/temp_p5/Run_2/data.nii',...
%    };
%
% slice timing parameters
%slice_params.slice_timing = false; % run slice timing correction
%slice_params.TR = 3; % repetition time in seconds
%slice_params.slice_order = 'descending'; % slice ordering (ascending or descending)
%
% fieldmap parameters
%field_params.fieldmap_undistortion = false; % run fieldmap undistortion
%field_params.fmap_magn = '';
%field_params.fmap_phase = '';
%field_params.fmap_te1 = 6.0; % shorter echo time in ms
%field_params.fmap_te2 = 7.02; % longer echo time in ms
%field_params.fmap_blipdir = -1; % phase-encoding direction
%field_params.fmap_BandwidthPerPixelPhaseEncode = 16.304; % phase-encoding BW Hz/px
%
% realignment parameters
%realign_params.unwarp = false;
%realign_params.mask = false;
%realign_params.c = [95 130 25];
%realign_params.r = [35 25 20];
%
% outlier parameters
%outlier_params.moco_out_mm_short = 0.4; % in mm
%outlier_params.moco_out_mm_long = 0.8; % in mm
%outlier_params.moco_out_deg_short = 0.5; % in deg
%outlier_params.moco_out_deg_long = 1.0; % in deg
%outlier_params.int_out_z = 2; % in z-score

% set spm default parameters
spm('defaults','FMRI');
spm_get_defaults('stats.maxmem',2^35); % maxmen indicates how much memory can be used
spm_get_defaults('cmdline',true); % no gui

% preprocessing summary
if length(img_input) > 1
    path_diagnosis = fullfile(fileparts(fileparts(img_input{1})),'diagnosis');
else
    path_diagnosis = fullfile(fileparts(img_input{1}),'diagnosis');
end

if ~exist(path_diagnosis,'dir') 
    mkdir(path_diagnosis); 
end

% create mask 
if realign_params.mask
    ft_cubic_mask(...
        img_input{1}, ...
        path_diagnosis, ...
        'moco_mask', ...
        realign_params.c, ...
        realign_params.r ...
        );
end

% path and filename of fieldmap phase
[path_fmap2, file_fmap2, ~] = fileparts(field_params.fmap_phase);

if slice_params.slice_timing
    for i = 1:length(img_input)
        
        % length of time series
        data_img = spm_vol(img_input{i});
        nt = length(data_img); % number of volumes
        nslices = data_img(1).dim(3); % number of slices
        
        scans = cell(nt,1);
        for j = 1:nt
            scans{j,1} = [img_input{i} ',' num2str(j)];
        end
        matlabbatch{1}.spm.temporal.st.scans = {scans};
        matlabbatch{1}.spm.temporal.st.nslices = nslices;
        matlabbatch{1}.spm.temporal.st.tr = slice_params.TR;
        matlabbatch{1}.spm.temporal.st.ta = slice_params.TR - slice_params.TR / nslices;
        if strcmp(slice_params.slice_order,'ascending')
            matlabbatch{1}.spm.temporal.st.so = 1:1:nslices;
            matlabbatch{1}.spm.temporal.st.refslice = 1;
        elseif strcmp(slice_params.slice_order,'descending')
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

if field_params.fieldmap_undistortion
       
    % total readout time
    total_readout = 1000/field_params.fmap_BandwidthPerPixelPhaseEncode; % in ms
     
    % calculate fieldmap
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = {field_params.fmap_phase};
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = {field_params.fmap_magn};
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [field_params.fmap_te1 field_params.fmap_te2];
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = field_params.fmap_blipdir;
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
cd(fileparts(img_input{1})); % change to first run directory
if realign_params.unwarp

    for i = 1:length(img_input)

        % get time series length
        data_img = spm_vol(img_input{i});
        nt = length(data_img);

        % get time series path and filename
        [path, file, ext] = fileparts(img_input{i});
        if slice_params.slice_timing
            file = ['a' file];
        end

        % get input data
        for j = 1:nt
            matlabbatch{1}.spm.spatial.realignunwarp.data(i).scans{j,1} = fullfile(path,[file ext ',' num2str(j)]);
            if field_params.fieldmap_undistortion
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

    % realignment (coregister only)
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 7;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
    if realign_params.mask
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = {fullfile(path_diagnosis,'moco_mask.nii')};
    else
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = {''};
    end

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

else

    for i = 1:length(img_input)

        % get time series length
        data_img = spm_vol(img_input{i});
        nt = length(data_img);

        % get time series path and filename
        [path, file, ext] = fileparts(img_input{i});
        if slice_params.slice_timing
            file = ['a' file];
        end

        % get input data
        for j = 1:nt
            matlabbatch{1}.spm.spatial.realign.estwrite.data{1,i}{j,1} = fullfile(path,[file ext ',' num2str(j)]);
        end

    end

    % realignment (coregister only)
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 7;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    if realign_params.mask
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {fullfile(path_diagnosis,'moco_mask.nii')};
    else
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {''};
    end
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 7;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'u';
    
end

% run
spm_jobman('run',matlabbatch);

% clear
clear matlabbatch

% check realignment processing
outlier_all = [];
for i = 1:length(img_input)

    % change to single run
    cd(fileparts(img_input{i}));
    
    [path, file, ext] = fileparts(img_input{i});
    if slice_params.slice_timing
        file = ['a' file];
    end
    
    % plot realignment parameters
    ft_plot_moco(...
        ['rp_' file '.txt'], ...
        'spm', ...
        path_diagnosis, ...
        ['rp_run_' num2str(i)] ...
        );
    
    % check for realignment and/or intensity outliers
    outlier = ft_outlier(...
        ['rp_' file '.txt'], ...
        ['u' file ext], ...
        outlier_params, ...
        'spm', ...
        fullfile(path,'outlier') ...
        );
    outlier_all = [outlier_all ; length(outlier)];
    
end

% open textfile
fileID = fopen(fullfile(path_diagnosis,['preprocessing_summary_u' file '.txt']),'w');
fprintf(fileID,'List of input parameters\n');
fprintf(fileID,'----------\n\n');
fprintf(fileID,'Preprocessed data\n');
for i = 1:length(img_input)
    fprintf(fileID,[img_input{i} '\n']);
end
fprintf(fileID,['run slice timing correction: ' mat2str(slice_params.slice_timing) '\n']);
fprintf(fileID,['slice timing parameter (TR): ' num2str(slice_params.TR) '\n']);
fprintf(fileID,['slice timing parameter (order): ' slice_params.slice_order '\n']);
fprintf(fileID,['run fieldmap undistortion: ' mat2str(field_params.fieldmap_undistortion) '\n']);
fprintf(fileID,['fieldmap undistortion parameter (magn): ' field_params.fmap_magn '\n']);
fprintf(fileID,['fieldmap undistortion parameter (phase): ' field_params.fmap_phase '\n']);
fprintf(fileID,['fieldmap undistortion parameter (TE1): ' num2str(field_params.fmap_te1) '\n']);
fprintf(fileID,['fieldmap undistortion parameter (TE2): ' num2str(field_params.fmap_te2) '\n']);
fprintf(fileID,['fieldmap undistortion parameter (blipdir): ' num2str(field_params.fmap_blipdir) '\n']);
fprintf(fileID,['fieldmap undistortion parameter (bw): ' num2str(field_params.fmap_BandwidthPerPixelPhaseEncode) '\n']);
fprintf(fileID,['run unwarp: ' mat2str(realign_params.unwarp) '\n']);
fprintf(fileID,['realign parameter (mask): ' mat2str(realign_params.mask) '\n']);
fprintf(fileID,['realign parameter (c): ' num2str(realign_params.c) '\n']);
fprintf(fileID,['realign parameter (r): ' num2str(realign_params.r) '\n']);
fprintf(fileID,'Percentage of within-run motion and intensity outliers\n');
fprintf(fileID,['motion threshold (mm, short): ' num2str(outlier_params.moco_out_mm_short) '\n']);
fprintf(fileID,['motion threshold (mm, long): ' num2str(outlier_params.moco_out_mm_long) '\n']);
fprintf(fileID,['motion threshold (deg, short): ' num2str(outlier_params.moco_out_deg_short) '\n']);
fprintf(fileID,['motion threshold (deg, long): ' num2str(outlier_params.moco_out_deg_long) '\n']);
fprintf(fileID,['intensity threshold (z-score): ' num2str(outlier_params.int_out_z) '\n']);
fprintf(fileID,'----------\n\n');

for i  = 1:length(img_input)

    % time series length
    data_img = spm_vol(img_input{i});
    nt = length(data_img);
    
    % get within-run outlier percentage
    outlier_percentage = outlier_all(i) / nt * 100;
    
    % get ratio of outliers within time series
    fprintf(fileID,'%.2f\n', outlier_percentage);
    
end

% close file
fclose(fileID);

% get time series of first volumes
data_img_out = spm_vol(img_input{1});
for i = 1:length(img_input)

    [path, file, ext] = fileparts(img_input{i});
    if slice_params.slice_timing
        file = ['a' file];
    end
    
    data_img = spm_vol(fullfile(path,['u' file ext]));
    data_array = spm_read_vols(data_img);
    
    data_img_out(i).fname = fullfile(path_diagnosis, ['vol1_u' file '.nii']);
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
    if slice_params.slice_timing
        file = ['a' file];
    end
    
    data_img = spm_vol(fullfile(path,['u' file ext]));
    data_array = spm_read_vols(data_img);
    
    data_mean_array = data_mean_array + sum(data_array,4);
    counter = counter + length(data_img);
    
end
data_mean_array = data_mean_array ./ counter;
data_img_out.fname = fullfile(path_diagnosis, ['mean_u' file '.nii']);
spm_write_vol(data_img_out, data_mean_array);
