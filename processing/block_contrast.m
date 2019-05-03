% GLM analysis

% This script computes a fixed-effects GLM for multiple runs within one
% session. The script is designed for a session with three different
% conditions, which are in my case two experimental ones and one baseline
% condition. Stimulus parameters are read from a spm12 compatible *.mat 
% file. No regressors of no interest are specified.

% created by Daniel Haenelt
% Date created: 10-12-2018
% Last modified: 03-05-2019

% input data
img_input = {
    '/nobackup/actinium2/haenelt/VasoTest/flicker/Run_1/ubold_basis_upsampled.nii',...
    '/nobackup/actinium2/haenelt/VasoTest/flicker/Run_2/ubold_basis_upsampled.nii',...
    '/nobackup/actinium2/haenelt/VasoTest/flicker/Run_3/ubold_basis_upsampled.nii',...
    '/nobackup/actinium2/haenelt/VasoTest/flicker/Run_4/ubold_basis_upsampled.nii',...
    '/nobackup/actinium2/haenelt/VasoTest/flicker/Run_5/ubold_basis_upsampled.nii',...
    '/nobackup/actinium2/haenelt/VasoTest/flicker/Run_6/ubold_basis_upsampled.nii',...
    };

cond_input = {
    '/nobackup/actinium2/haenelt/VasoTest/flicker/Run_1/logfiles/VasoTest_flicker_Run1_Cond.mat',...
    '/nobackup/actinium2/haenelt/VasoTest/flicker/Run_2/logfiles/VasoTest_flicker_Run2_Cond.mat',...
    '/nobackup/actinium2/haenelt/VasoTest/flicker/Run_3/logfiles/VasoTest_flicker_Run3_Cond.mat',...
    '/nobackup/actinium2/haenelt/VasoTest/flicker/Run_4/logfiles/VasoTest_flicker_Run4_Cond.mat',...
    '/nobackup/actinium2/haenelt/VasoTest/flicker/Run_5/logfiles/VasoTest_flicker_Run5_Cond.mat',...
    '/nobackup/actinium2/haenelt/VasoTest/flicker/Run_6/logfiles/VasoTest_flicker_Run6_Cond.mat',...
    };

% parameters
TR = 2.5; % repetition time  in s
cutoff_highpass = 100; % 1/cutoff_highpass frequency in Hz
nconds = 2; % only 2 and 3 are supported
name_output = 'bold'; % basename of output contrasts
output_folder = 'contrast_bold'; % name of folder where spm.mat is saved

% add spm to path
pathSPM = '/data/pt_01880/source/spm12'; 

%%% do not edit below %%%

% add paths to the interpreter's search path
addpath(pathSPM);
spm('defaults','FMRI');
spm_get_defaults('stats.maxmem',2^35); % maxmen indicates how much memory can be used
spm_get_defaults('cmdline',true); % no gui

% output folder is taken from the first entry of the input list
path_output = fullfile(fileparts(fileparts(img_input{1})),output_folder);
if ~exist(path_output,'dir')
    mkdir(path_output);
end
cd(path_output);

% session name is taken from the first entry of the input list
[~, name_sess,~] = fileparts(fileparts(fileparts(img_input{1})));

% number of volumes is taken from the first entry of the input list
data_img = spm_vol(img_input{1});
nt = length(data_img);

% fmri model
matlabbatch{1}.spm.stats.fmri_spec.dir = {path_output};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

for i = 1:length(img_input)    
    for j = 1:nt
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans{j,1} = [img_input{i} ',' num2str(j)];
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi = {cond_input{i}};
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf = cutoff_highpass;    
end

% hrf model
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.3;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

% run
spm_jobman('run',matlabbatch);

% clear matlabbatch
clear matlabbatch

% model estimation
matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(path_output,'SPM.mat')};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

% run
spm_jobman('run',matlabbatch);

% clear after completion
clear matlabbatch

% calculate contrast
% condition names are taken from the first entry of the condition file list
if nconds == 2
    
    cond = load(cond_input{1});
    name_contrast = {
        [cond.names{2} '_' cond.names{1}],...
        [cond.names{1} '_' cond.names{2}],...
        };

    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(path_output,'SPM.mat')};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = name_contrast{1};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [-1 1];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = name_contrast{2};
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [1 -1];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'replsc';
    matlabbatch{1}.spm.stats.con.delete = 0; % append contrast to old contrasts
    
elseif nconds == 3
    
    cond = load(cond_input{1});
    name_contrast = {
        [cond.names{3} '_' cond.names{2}],...
        [cond.names{2} '_' cond.names{3}],...
        [cond.names{3} '_' cond.names{1}],...
        [cond.names{2} '_' cond.names{1}],...
        };

    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(path_output,'SPM.mat')};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = name_contrast{1};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [0 -1 1];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = name_contrast{2};
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [0 1 -1];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'replsc';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = name_contrast{3};
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [-1 0 1];
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'replsc';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = name_contrast{4};
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [-1 1 0];
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'replsc';
    matlabbatch{1}.spm.stats.con.delete = 0; % append contrast to old contrasts

end

% run
spm_jobman('run',matlabbatch);

% clear after completion
clear matlabbatch

% spmT images
path_in = path_output;
path_out = fullfile(fileparts(fileparts(fileparts(img_input{1}))),'results','spmT','native'); 
if ~exist(path_out,'dir')
    mkdir(path_out);
end

% sort results
for i = 1:length(name_contrast)
    % copy file
    copyfile(...
        fullfile(path_in,['spmT_' sprintf('%04d',i) '.nii']),...
        fullfile(path_out,['spmT_' name_output '_' name_contrast{i} '_' name_sess '.nii']));
    
    % set nan to zero
    data_img = spm_vol(fullfile(path_out,['spmT_' name_output '_' name_contrast{i} '_' name_sess '.nii']));
    data_array = spm_read_vols(data_img);
    data_array(isnan(data_array)) = 0;
    spm_write_vol(data_img,data_array);
end    

% con images
path_in = path_output;
path_out = fullfile(fileparts(fileparts(fileparts(img_input{1}))),'results','con','native'); 
if ~exist(path_out,'dir')
    mkdir(path_out);
end

% sort results
for i = 1:length(name_contrast)
    % copy file
    copyfile(...
        fullfile(path_in,['con_' sprintf('%04d',i) '.nii']),...
        fullfile(path_out,['con_' name_output '_' name_contrast{i} '_' name_sess '.nii']));
    
    % set nan to zero
    data_img = spm_vol(fullfile(path_out,['con_' name_output '_' name_contrast{i} '_' name_sess '.nii']));
    data_array = spm_read_vols(data_img);
    data_array(isnan(data_array)) = 0;
    spm_write_vol(data_img,data_array);
end