% GLM analysis

% This script computes a fixed-effects GLM for one or multiple runs within 
% one session. The script is designed for paradigms with 2-4 experimental 
% conditions. However, for more than 2 conditions, not all contrasts are 
% computed:) Stimulus parameters are read from a spm12 compatible *.mat 
% file. Optionally, regressors of no interest can be specified.

% created by Daniel Haenelt
% Date created: 10-12-2018
% Last modified: 05-09-2019

% input data
img_input = {
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_1/udata.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_2/udata.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_3/udata.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_4/udata.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_5/udata.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_6/udata.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_7/udata.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_8/udata.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_9/udata.nii',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_10/udata.nii',...
    };

cond_input = {
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_1/logfiles/p3_GE_EPI1_Run1_nonrivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_2/logfiles/p3_GE_EPI1_Run2_nonrivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_3/logfiles/p3_GE_EPI1_Run3_nonrivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_4/logfiles/p3_GE_EPI1_Run4_nonrivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_5/logfiles/p3_GE_EPI1_Run5_nonrivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_6/logfiles/p3_GE_EPI1_Run6_nonrivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_7/logfiles/p3_GE_EPI1_Run7_nonrivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_8/logfiles/p3_GE_EPI1_Run8_nonrivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_9/logfiles/p3_GE_EPI1_Run9_nonrivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_10/logfiles/p3_GE_EPI1_Run10_nonrivalry_Cond.mat',...
    };

multi_input = {
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_1/logfiles/correlation_regressor_Run1.txt',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_2/logfiles/correlation_regressor_Run2.txt',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_3/logfiles/correlation_regressor_Run3.txt',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_4/logfiles/correlation_regressor_Run4.txt',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_5/logfiles/correlation_regressor_Run5.txt',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_6/logfiles/correlation_regressor_Run6.txt',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_7/logfiles/correlation_regressor_Run7.txt',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_8/logfiles/correlation_regressor_Run8.txt',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_9/logfiles/correlation_regressor_Run9.txt',...
    '/data/pt_01880/Experiment2_Rivalry/p3/nonrivalry/GE_EPI1/Run_10/logfiles/correlation_regressor_Run10.txt',...
    };

% parameters
TR = 3; % repetition time  in s
cutoff_highpass = 180; % 1/cutoff_highpass frequency in Hz
nconds = 3; % only 2-4 are supported
name_sess = 'test3'; % name of session (if multiple sessions exist)
name_output = ''; % basename of output contrasts
output_folder = 'contrast3'; % name of folder where spm.mat is saved
run_multi = true;

% lowpass filter of time series
lowpass = false;
cutoff_lowpass = 0;
order_lowpass = 3;

% add spm and lib to path
pathSPM = '/data/pt_01880/source/spm12'; 
pathLIB = '/home/raid2/haenelt/projects/scripts/lib/';

%%% do not edit below %%%

% add paths to the interpreter's search path
addpath(pathSPM);
spm('defaults','FMRI');
spm_get_defaults('stats.maxmem',2^35); % maxmen indicates how much memory can be used
spm_get_defaults('cmdline',true); % no gui
addpath(fullfile(pathLIB,'processing'));

% lowpass filtering
if lowpass
    for i = 1:length(img_input)
        lowpass_filter(img_input{i}, TR, cutoff_lowpass, order_lowpass, pathSPM);
    end
    
    % change input to lowpass filtered time series
    [filepath, name, ext] = fileparts(img_input{i});
    img_input{i} = fullfile(filepath,['l' name ext]);
end

% output folder is taken from the first entry of the input list
if length(img_input) == 1
    path_output = fullfile(fileparts(img_input{1}),output_folder);
else
    path_output = fullfile(fileparts(fileparts(img_input{1})),output_folder);
end
if ~exist(path_output,'dir')
    mkdir(path_output);
end
cd(path_output);

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
    if run_multi
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = {multi_input{i}};
    else
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = {''};
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf = cutoff_highpass;    
end

% hrf model
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.3; % -Inf for no implicit mask
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

elseif nconds == 4
    
    cond = load(cond_input{1});
    name_contrast = {
        [cond.names{1} '_' cond.names{2}],...
        [cond.names{2} '_' cond.names{1}],...
        [cond.names{1} '_' cond.names{3}],...
        [cond.names{2} '_' cond.names{4}],...
        };

    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(path_output,'SPM.mat')};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = name_contrast{1};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1 0 0];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = name_contrast{2};
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [-1 1 0 0];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'replsc';
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = name_contrast{3};
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [1 0 -1 0];
    matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'replsc';
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = name_contrast{4};
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [0 1 0 -1];
    matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'replsc';
    matlabbatch{1}.spm.stats.con.delete = 0; % append contrast to old contrasts
    
end

% run
spm_jobman('run',matlabbatch);

% clear after completion
clear matlabbatch

% make folders for spmT and con images
if length(img_input) == 1
    path_spmT = fullfile(fileparts(img_input{1}),'results','spmT','native'); 
    path_con = fullfile(fileparts(img_input{1}),'results','con','native'); 
else
    path_spmT = fullfile(fileparts(fileparts(fileparts(img_input{1}))),'results','spmT','native'); 
    path_con = fullfile(fileparts(fileparts(fileparts(img_input{1}))),'results','con','native'); 
end

if ~exist(path_spmT,'dir')
    mkdir(path_spmT);
end
if ~exist(path_con,'dir')
    mkdir(path_con);
end

% sort results
for i = 1:length(name_contrast)
    % copy file
    file_in_spmT = fullfile(path_output,['spmT_' sprintf('%04d',i) '.nii']);
    file_in_con = fullfile(path_output,['con_' sprintf('%04d',i) '.nii']);
    
    if isempty(name_output) && isempty(name_sess)
        file_out_spmT = fullfile(path_spmT,['spmT_' name_contrast{i} '.nii']);
        file_out_con = fullfile(path_con,['con_' name_contrast{i} '.nii']);
    elseif isempty(name_output) && ~isempty(name_sess)
        file_out_spmT = fullfile(path_spmT,['spmT_' name_contrast{i} '_' name_sess '.nii']);
        file_out_con = fullfile(path_con,['con_' name_contrast{i} '_' name_sess '.nii']);
    elseif ~isempty(name_output) && isempty(name_sess)
        file_out_spmT = fullfile(path_spmT,['spmT_' name_output '_' name_contrast{i} '.nii']);
        file_out_con = fullfile(path_con,['con_' name_output '_' name_contrast{i} '.nii']);  
    else
        file_out_spmT = fullfile(path_spmT,['spmT_' name_output '_' name_contrast{i} '_' name_sess '.nii']);
        file_out_con = fullfile(path_con,['con_' name_output '_' name_contrast{i} '_' name_sess '.nii']);        
    end
    copyfile(file_in_spmT, file_out_spmT);
    copyfile(file_in_con, file_out_con);
    
    % set nan to zero
    data_img = spm_vol(file_out_spmT);
    data_array = spm_read_vols(data_img);
    data_array(isnan(data_array)) = 0;
    spm_write_vol(data_img,data_array);
    
    data_img = spm_vol(file_out_con);
    data_array = spm_read_vols(data_img);
    data_array(isnan(data_array)) = 0;
    spm_write_vol(data_img,data_array);
end
