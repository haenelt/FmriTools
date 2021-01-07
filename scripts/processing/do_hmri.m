% hMRI processing

% Computation of quantitative maps from acquired MPMs. The script is
% based on the batch GUI. All input files are defined as cell arrays. The 
% frame number in the filename string is omitted since 1 (3d volume) is the 
% default parameter for file selections. The script consists of the 
% following steps: (1) 7T configurations are loaded. (2) Images are 
% oriented to align with MNI space. Input for registration is the first pd
% echo. The registration is then applied to all other images (pdw, t1w, b1,
% b0). (3) Quantitative maps are created. (4) Resulting maps are reoriented
% back to native space. (5) Maps are scaled.

% created by Daniel Haenelt
% Date created: 11-09-2020
% Last modified: 07-01-2021

% input data
file_t1w = {
    '/data/pt_01880/test_mpm/t1_kp_mtflash3d_v1ax_0p5_0007/s2803151-123034-00001-00352-1.nii'
    '/data/pt_01880/test_mpm/t1_kp_mtflash3d_v1ax_0p5_0007/s2803151-123034-00001-00704-2.nii'
    '/data/pt_01880/test_mpm/t1_kp_mtflash3d_v1ax_0p5_0007/s2803151-123034-00001-01056-3.nii'
    '/data/pt_01880/test_mpm/t1_kp_mtflash3d_v1ax_0p5_0007/s2803151-123034-00001-01408-4.nii'
    '/data/pt_01880/test_mpm/t1_kp_mtflash3d_v1ax_0p5_0007/s2803151-123034-00001-01760-5.nii'
    '/data/pt_01880/test_mpm/t1_kp_mtflash3d_v1ax_0p5_0007/s2803151-123034-00001-02112-6.nii'
    };

file_pdw = {
    '/data/pt_01880/test_mpm/pd_kp_mtflash3d_v1ax_0p5_0008/s2803151-124913-00001-00352-1.nii'
    '/data/pt_01880/test_mpm/pd_kp_mtflash3d_v1ax_0p5_0008/s2803151-124913-00001-00704-2.nii'
    '/data/pt_01880/test_mpm/pd_kp_mtflash3d_v1ax_0p5_0008/s2803151-124913-00001-01056-3.nii'
    '/data/pt_01880/test_mpm/pd_kp_mtflash3d_v1ax_0p5_0008/s2803151-124913-00001-01408-4.nii'
    '/data/pt_01880/test_mpm/pd_kp_mtflash3d_v1ax_0p5_0008/s2803151-124913-00001-01760-5.nii'
    '/data/pt_01880/test_mpm/pd_kp_mtflash3d_v1ax_0p5_0008/s2803151-124913-00001-02112-6.nii'
    };

file_b0 = {
    '/data/pt_01880/test_mpm/gre_field_mapping_2mm_0005/s2803151-122403-00001-00001-1.nii'
    '/data/pt_01880/test_mpm/gre_field_mapping_2mm_0005/s2803151-122404-00001-00001-2.nii'
    '/data/pt_01880/test_mpm/gre_field_mapping_2mm_0006/s2803151-122404-00001-00001-2.nii'
    };

file_b1 = {
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122057-00001-00001-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122057-00001-00049-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122109-00002-00097-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122109-00002-00145-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122121-00003-00193-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122121-00003-00241-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122133-00004-00289-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122133-00004-00337-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122145-00005-00385-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122145-00005-00433-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122157-00006-00481-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122157-00006-00529-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122209-00007-00577-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122209-00007-00625-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122221-00008-00673-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122221-00008-00721-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122233-00009-00769-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122234-00009-00817-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122246-00010-00865-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122246-00010-00913-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122258-00011-00961-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122258-00011-01009-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122310-00012-01057-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122310-00012-01105-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122322-00013-01153-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122322-00013-01201-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122334-00014-01249-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122334-00014-01297-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122346-00015-01345-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2803151-122346-00015-01393-2.nii'
    };

% path to spm12 toolbox and hMRI configuration files.
pathHMRI_defaults = '/data/p_gr_weiskopf_software/spm12/toolbox/hMRI-cbs/config/local/hmri_CBS_7T_defaults.m';
pathHMRI_b1_defaults = '/data/p_gr_weiskopf_software/spm12/toolbox/hMRI-cbs/config/local/hmri_b1_CBS_7T_defaults.m';

%%% do not edit below %%%

% get SPM12 root directory from matlab search path
pathSPM = what('spm12').path;

% change current working directory. The working directoy is the parent
% folder of the pdw images.
path_output = fileparts(file_pdw{1});
cd(path_output);

% run toolbox configuration
matlabbatch{1}.spm.tools.hmri.hmri_config.hmri_setdef.customised = {pathHMRI_defaults};

% run auto-reorient
file_other = {file_pdw(2:end), file_t1w, file_b1, file_b0};
file_other = cat(1, file_other{:});
matlabbatch{2}.spm.tools.hmri.autoreor.reference = file_pdw(1);
matlabbatch{2}.spm.tools.hmri.autoreor.template = {fullfile(pathSPM, 'canonical/avg152T1.nii')};
matlabbatch{2}.spm.tools.hmri.autoreor.other = file_other;
matlabbatch{2}.spm.tools.hmri.autoreor.output.indir = 'yes';
matlabbatch{2}.spm.tools.hmri.autoreor.dep = 'individual';

% run create maps
matlabbatch{3}.spm.tools.hmri.create_mpm.subj.output.indir = 'yes';
matlabbatch{3}.spm.tools.hmri.create_mpm.subj.sensitivity.RF_us = '-';
matlabbatch{3}.spm.tools.hmri.create_mpm.subj.b1_type.i3D_EPI.b1input = file_b1;
matlabbatch{3}.spm.tools.hmri.create_mpm.subj.b1_type.i3D_EPI.b0input = file_b0;
matlabbatch{3}.spm.tools.hmri.create_mpm.subj.b1_type.i3D_EPI.b1parameters.b1defaults = {pathHMRI_b1_defaults};
matlabbatch{3}.spm.tools.hmri.create_mpm.subj.raw_mpm.MT = '';
matlabbatch{3}.spm.tools.hmri.create_mpm.subj.raw_mpm.PD = file_pdw;
matlabbatch{3}.spm.tools.hmri.create_mpm.subj.raw_mpm.T1 = file_t1w;
matlabbatch{3}.spm.tools.hmri.create_mpm.subj.popup = false;

% start job
spm_jobman('run', matlabbatch);

clear matlabbatch

% get transformation matrix from json
fid = fopen(fullfile(path_output, 'AutoReorient_output.json')); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid);
invM = jsondecode(str).invM;

% reorient
[~, basename, ext] = fileparts(file_pdw{1});
file_result = {
    fullfile(path_output, 'Results', [basename '_PD' ext])
    fullfile(path_output, 'Results', [basename '_R1' ext])
    fullfile(path_output, 'Results', [basename '_R2s_WOLS' ext])
    };
matlabbatch{1}.spm.util.reorient.srcfiles = file_result;
matlabbatch{1}.spm.util.reorient.transform.transM = invM;
matlabbatch{1}.spm.util.reorient.prefix = 'r';

% start job
spm_jobman('run', matlabbatch);

clear matlabbatch

% scale and invert results (pd, r1, r2s)
file_result = {
    fullfile(path_output, 'Results', ['r' basename '_PD' ext])
    fullfile(path_output, 'Results', ['r' basename '_R1' ext])
    fullfile(path_output, 'Results', ['r' basename '_R2s_WOLS' ext])
    };

for i = 1:length(file_result)
    data_img = spm_vol(file_result{i});
    data_array = spm_read_vols(data_img);
    
    if strfind(file_result{i},'PD')
        val_min = 0;
        val_max = 100;
        val_invert = true;
    elseif strfind(file_result{i},'R1')
        val_min = 0;
        val_max = 5;
        val_invert = false;
    elseif strfind(file_result{i}, 'R2s')
        val_min = 0;
        val_max = 250;
        val_invert = false;
    end
    
    data_array(data_array < val_min) = val_min;
    data_array(data_array > val_max) = val_max;
    
    if val_invert
        data_array = max(data_array) - data_array;
    end
    
    [path, basename, ext] = fileparts(file_result{i});
    if val_invert
        data_img.fname = fullfile(path, [basename '_scaled_inverted' ext]);
    else
        data_img.fname = fullfile(path, [basename '_scaled' ext]);
    end
    spm_write_vol(data_img, data_array);
end
