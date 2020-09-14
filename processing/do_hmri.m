% hMRI processing

% Computation of quantitative maps from acquired MPMs. The script is
% completely based on the batch GUI. All input files are defined as cell
% arrays. The frame number in the filename string is omitted since 1 (3d 
% volume) is the default parameter for file selections. The script consists 
% of the following steps: (1) 7T configurations are loaded. (2) Images are
% oriented to align with MNI space. Input for registration is the first pd
% echo. The registration is then applied to all other images (pdw, t1w, b1,
% b0). (3) Quantitative maps are created.

% created by Daniel Haenelt
% Date created: 11-09-2020
% Last modified: 14-09-2020

% input data
file_t1w = {
    '/data/pt_01880/test_mpm/t1_kp_mtflash3d_v1ax_0p5_0007/s2899875-122853-00001-00352-1.nii'
    '/data/pt_01880/test_mpm/t1_kp_mtflash3d_v1ax_0p5_0007/s2899875-122853-00001-00704-2.nii'
    '/data/pt_01880/test_mpm/t1_kp_mtflash3d_v1ax_0p5_0007/s2899875-122853-00001-01056-3.nii'
    '/data/pt_01880/test_mpm/t1_kp_mtflash3d_v1ax_0p5_0007/s2899875-122853-00001-01408-4.nii'
    '/data/pt_01880/test_mpm/t1_kp_mtflash3d_v1ax_0p5_0007/s2899875-122853-00001-01760-5.nii'
    '/data/pt_01880/test_mpm/t1_kp_mtflash3d_v1ax_0p5_0007/s2899875-122853-00001-02112-6.nii'
    };

file_pdw = {
    '/data/pt_01880/test_mpm/pd_kp_mtflash3d_v1ax_0p5_0008/s2899875-124739-00001-00352-1.nii'
    '/data/pt_01880/test_mpm/pd_kp_mtflash3d_v1ax_0p5_0008/s2899875-124739-00001-00704-2.nii'
    '/data/pt_01880/test_mpm/pd_kp_mtflash3d_v1ax_0p5_0008/s2899875-124739-00001-01056-3.nii'
    '/data/pt_01880/test_mpm/pd_kp_mtflash3d_v1ax_0p5_0008/s2899875-124739-00001-01408-4.nii'
    '/data/pt_01880/test_mpm/pd_kp_mtflash3d_v1ax_0p5_0008/s2899875-124739-00001-01760-5.nii'
    '/data/pt_01880/test_mpm/pd_kp_mtflash3d_v1ax_0p5_0008/s2899875-124739-00001-02112-6.nii'
    };

file_b0 = {
    '/data/pt_01880/test_mpm/gre_field_mapping_2mm_0005/s2899875-122129-00001-00001-1.nii'
    '/data/pt_01880/test_mpm/gre_field_mapping_2mm_0005/s2899875-122130-00001-00001-2.nii'
    '/data/pt_01880/test_mpm/gre_field_mapping_2mm_0006/s2899875-122130-00001-00001-2.nii'
    };

file_b1 = {
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-121823-00001-00001-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-121823-00001-00049-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-121835-00002-00097-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-121835-00002-00145-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-121848-00003-00193-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-121848-00003-00241-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-121900-00004-00289-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-121900-00004-00337-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-121912-00005-00385-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-121912-00005-00433-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-121924-00006-00481-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-121924-00006-00529-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-121936-00007-00577-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-121936-00007-00625-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-121948-00008-00673-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-121948-00008-00721-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-122000-00009-00769-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-122000-00009-00817-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-122012-00010-00865-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-122012-00010-00913-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-122024-00011-00961-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-122024-00011-01009-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-122036-00012-01057-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-122036-00012-01105-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-122048-00013-01153-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-122048-00013-01201-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-122100-00014-01249-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-122100-00014-01297-2.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-122112-00015-01345-1.nii'
    '/data/pt_01880/test_mpm/kp_seste_b1map_v1a_TM34910_0004/s2899875-122112-00015-01393-2.nii'
    };

% path to spm12 toolbox and hMRI configuration files.
pathSPM = '/data/p_gr_weiskopf_software/spm12';
pathHMRI_defaults = '/data/pt_01880/temp_scripts/stripes/p1/hmri/hmri_CBS_7T_defaults_daniel.m';
pathHMRI_b1_defaults = '/data/pt_01880/temp_scripts/stripes/p1/hmri/hmri_b1_CBS_7T_defaults_daniel.m';

%%% do not edit below %%%

% add spm to matlab search path
addpath(pathSPM);

% change current working directory. The working directoy is the parent
% folder of the pdw images.
cd(fileparts(fileparts(file_pdw{1})));

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
matlabbatch{3}.spm.tools.hmri.create_mpm.subj.popup = true;

% start job
spm_jobman('run', matlabbatch);