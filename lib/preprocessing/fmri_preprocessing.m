function fmri_preprocessing(img_input, slice_params, field_params, pathSPM)
% This function performs slice time correction, fieldmap undistortion and
% motion correction in the SPM12 framework which can be applied to a
% session consisting of multiple runs. Slice time correction and fieldmap 
% undistortion are optional.
% Fieldmap undistortion is performed without skullstripping. The total
% readout time is taken as the inverse of the BandwidthPerPixelPhaseEncode
% in Hz/px which can be found in the dicom tag (0019, 1028).
% 
% Inputs:
    % input: cell array of filenames of input time series.
    % slice_params: struct of slice time correction parameters.
    % field_params: struct of fieldmap parameters.
    % pathSPM: path to spm12 folder.

% created by Daniel Haenelt
% Date created: 26-02-2019
% Last modified: 19-02-2020

% add spm to path
addpath(pathSPM);
spm('defaults','FMRI');
spm_get_defaults('stats.maxmem',2^35);% maxmen indicates how much memory can be used
spm_get_defaults('cmdline',true);% no gui

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