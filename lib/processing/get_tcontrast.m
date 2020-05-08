function get_tcontrast(cond_input, path_contrast, name_output, name_sess, hrf_derivative, ...
    regressor_nointerest, pathSPM)

% This function calculates contrasts (con, spmT) and percent signal changes
% (psc) from a general linear model using spm. Only 2-4 experimental
% conditions are supported in a block design. Not all contrast permutations
% are computed. Percent signal change calculations are largely based on
% Mazaika 2009. Optionally, contrast estimates are computed with
% consideration of the hrf amplitude and its temporal derivative in the
% design matrix. To mitigate latency-induced amplitude bias in the
% resulting effect sizes, percent signal changes are also computed based on
% the method by Calhoun et al. 2004.

% Inputs:
    % cond_input: cell array containing filenames of condition files.
    % path_contrast: path where spm.mat is saved.
    % name_output: name of output (optional).
    % name_sess: name of session (optional).
    % hrf_derivative: use hrf temporal derivatives (boolean).
    % regressor_nointerest: include number of regressors of no interest.
    % pathSPM: path to spm12 folder.
    
% created by Daniel Haenelt
% Date created: 13-03-2020
% Last modified: 08-05-2020

if ~exist('name_output', 'var')  
    name_output = '';
end

if ~exist('name_sess', 'var')  
    name_sess = '';
end

if ~exist('hrf_derivative', 'var')  
    hrf_derivative = false;
end

% add spm to path
addpath(pathSPM);

% print to console
disp('compute percent signal changes');

% check if noise regressor vector has correct length
if length(cond_input) ~= length(regressor_nointerest)
    error('regressor_nointerest has wrong size!');
end

% load conditions first cell entry
cond = load(cond_input{1});

% get number of conditions and number of runs
nruns = length(cond_input);
nconds = length(cond.names);

% make output folders
if ~isempty(name_sess)
    path_spmT = fullfile(fileparts(fileparts(path_contrast)),'results','spmT','native');
    path_con = fullfile(fileparts(fileparts(path_contrast)),'results','con','native');
    path_psc = fullfile(fileparts(fileparts(path_contrast)),'results','psc','native');
else
    path_spmT = fullfile(fileparts(path_contrast),'results','spmT','native'); 
    path_con = fullfile(fileparts(path_contrast),'results','con','native'); 
    path_psc = fullfile(fileparts(path_contrast),'results','psc','native');
end

if ~exist(path_spmT,'dir')
    mkdir(path_spmT);
end
if ~exist(path_con,'dir')
    mkdir(path_con);
end
if ~exist(path_psc,'dir')
    mkdir(path_psc);
end

% get condition names
if nconds == 2
    name_contrast = {
        [cond.names{2} '_' cond.names{1}],...
        [cond.names{1} '_' cond.names{2}],...
        cond.names{1},...
        cond.names{2},...
        };
elseif nconds == 3
    name_contrast = {
        [cond.names{3} '_' cond.names{2}],...
        [cond.names{2} '_' cond.names{3}],...
        [cond.names{3} '_' cond.names{1}],...
        [cond.names{2} '_' cond.names{1}],...
        cond.names{1},...
        cond.names{2},...
        cond.names{3},...
        };
elseif nconds == 4
    name_contrast = {
        [cond.names{1} '_' cond.names{2}],...
        [cond.names{2} '_' cond.names{1}],...
        [cond.names{1} '_' cond.names{3}],...
        [cond.names{2} '_' cond.names{4}],...
        cond.names{1},...
        cond.names{2},...
        cond.names{3},...
        cond.names{4},...
        };
end

% double contrast names if design matrix with temporal derivative
if hrf_derivative == true
    name_contrast = repelem(name_contrast,2);
    name_contrast(2:2:end) = cellfun(@(c)[c '_diff'], name_contrast(2:2:end),'uni',false);
end

% get basenames for output files
basename = {};
for i = 1:length(name_contrast)
    if isempty(name_output) && isempty(name_sess)
        basename{i} = name_contrast{i};
    elseif isempty(name_output) && ~isempty(name_sess)
        basename{i} =  [name_contrast{i} '_' name_sess];
    elseif ~isempty(name_output) && isempty(name_sess)
        basename{i} = [name_output '_' name_contrast{i}];
    else
        basename{i} = [name_output '_' name_contrast{i} '_' name_sess];
    end
end

% get contrast vector
if hrf_derivative == true
    if nconds == 2
        c = [-1 0 1 0 ; 
             0 -1 0 1 ; 
             1 0 -1 0 ; 
             0 1 0 -1 ; 
             1 0 0 0 ; 
             0 1 0 0 ; 
             0 0 1 0 ; 
             0 0 0 1];
    elseif nconds == 3
        c = [0 0 -1 0 1 0 ; 
             0 0 0 -1 0 1 ; 
             0 0 1 0 -1 0 ; 
             0 0 0 1 0 -1 ; 
             -1 0 0 0 1 0 ; 
             0 -1 0 0 0 1 ; 
             -1 0 1 0 0 0 ; 
             0 -1 0 1 0 0 ; 
             1 0 0 0 0 0 ; 
             0 1 0 0 0 0 ; 
             0 0 1 0 0 0 ; 
             0 0 0 1 0 0 ; 
             0 0 0 0 1 0 ;
             0 0 0 0 0 1];
    elseif nconds == 4
        c = [1 0 -1 0 0 0 0 0 ; 
             0 1 0 -1 0 0 0 0 ; 
             -1 0 1 0 0 0 0 0 ; 
             0 -1 0 1 0 0 0 0 ; 
             1 0 0 0 -1 0 0 0 ; 
             0 1 0 0 0 -1 0 0 ; 
             0 0 1 0 0 0 -1 0 ; 
             0 0 0 1 0 0 0 -1 ; 
             1 0 0 0 0 0 0 0 ; 
             0 1 0 0 0 0 0 0 ; 
             0 0 1 0 0 0 0 0 ; 
             0 0 0 1 0 0 0 0 ; 
             0 0 0 0 1 0 0 0 ; 
             0 0 0 0 0 1 0 0 ; 
             0 0 0 0 0 0 1 0 ; 
             0 0 0 0 0 0 0 1];
    end
else
    if nconds == 2
        c = [-1 1 ; 
             1 -1 ; 
             1 0 ; 
             0 1];
    elseif nconds == 3
        c = [0 -1 1 ; 
             0 1 -1 ; 
             -1 0 1 ; 
             -1 1 0 ; 
             1 0 0 ; 
             0 1 0 ; 
             0 0 1];
    elseif nconds == 4
        c = [1 -1 0 0 ; 
             -1 1 0 0 ; 
             1 0 -1 0 ; 
             0 1 0 -1 ; 
             1 0 0 0 ; 
             0 1 0 0 ; 
             0 0 1 0 ; 
             0 0 0 1];
    end
end

% add for multiple runs and scale
contrast = [];
for i=1:nruns
    contrast = [contrast c/nruns zeros(length(c),regressor_nointerest(i))];
end

% add zeros for constant term
contrast = [contrast zeros(size(contrast,1),nruns)];

% compute tscores
matlabbatch{1}.spm.stats.con.spmmat = {fullfile(path_contrast,'SPM.mat')};
matlabbatch{1}.spm.stats.con.delete = 0; % append contrast to old contrasts
for i = 1:length(name_contrast)
    matlabbatch{1}.spm.stats.con.consess{i}.tcon.name = name_contrast{i};
    matlabbatch{1}.spm.stats.con.consess{i}.tcon.weights = contrast(i,:);
    matlabbatch{1}.spm.stats.con.consess{i}.tcon.sessrep = 'none';
end

% run
spm_jobman('run',matlabbatch);

% clear after completion
clear matlabbatch

% sort results to output folder
for i = 1:length(name_contrast)

    % copy file
    file_in_spmT = fullfile(path_contrast,['spmT_' sprintf('%04d',i) '.nii']);
    file_in_con = fullfile(path_contrast,['con_' sprintf('%04d',i) '.nii']);
    
    file_out_spmT = fullfile(path_spmT,['spmT_' basename{i} '.nii']);
    file_out_con = fullfile(path_con,['con_' basename{i} '.nii']);
    
    copyfile(file_in_spmT, file_out_spmT);
    copyfile(file_in_con, file_out_con);
    
    % set nan to zero
    data_img = spm_vol(file_out_spmT);
    data_array = spm_read_vols(data_img);
    data_array(isnan(data_array)) = 0;
    spm_write_vol(data_img, data_array);
    
    data_img = spm_vol(file_out_con);
    data_array = spm_read_vols(data_img);
    data_array(isnan(data_array)) = 0;
    spm_write_vol(data_img, data_array);
    
end

% compute psc

% load mask
mask_img = spm_vol(fullfile(path_contrast,'mask.nii'));
mask_array = spm_read_vols(mask_img);
mask_array(isnan(mask_array)) = 0;

% load constant
nbeta = size(contrast,2);
constant_array = zeros(size(mask_array));
for i = 1:nruns
    file_constant = fullfile(path_contrast,['beta_' sprintf('%04d',nbeta-nruns+i) '.nii']);
    constant_img = spm_vol(file_constant);
    constant_array = constant_array + spm_read_vols(constant_img);
end
constant_array(isnan(constant_array)) = 0;
constant_array = constant_array / nruns;

% get baseline value
bmean = mean(constant_array(mask_array > 0));

% get peak value
load(fullfile(path_contrast,'SPM.mat'))
peak = max(SPM.xX.X(:));

% print to console
disp(['Peak value for psc calculation: ' num2str(peak)]);
disp(['Baseline value for psc calculation: ' num2str(bmean)]);

for i = 1:length(name_contrast)
    
    file_in_con = fullfile(path_con,['con_' basename{i} '.nii']);
    file_out_psc = fullfile(path_psc,['psc_' basename{i} '.nii']);
    
    psc_img = spm_vol(file_in_con);
    psc_array = spm_read_vols(psc_img);
    psc_array = psc_array * peak / bmean * 100;
    psc_img.fname = file_out_psc;
    spm_write_vol(psc_img, psc_array);
    
end

% psc free of latency-induced amplitude bias
if hrf_derivative == true
    for i = 1:2:length(name_contrast)

        file_in_psc1 = fullfile(path_psc,['psc_' basename{i} '.nii']);
        file_in_psc2 = fullfile(path_psc,['psc_' basename{i+1} '.nii']);
        file_out_psc = fullfile(path_psc,['psc_' basename{i} '_calhoun.nii']);

        beta_img = spm_vol(file_in_psc1);
        beta_diff_img = spm_vol(file_in_psc2);
        beta_array = spm_read_vols(beta_img);
        beta_diff_array = spm_read_vols(beta_diff_img);

        beta_array = sign(beta_array) .* sqrt(beta_array.^2 + beta_diff_array.^2);
        beta_img.fname = file_out_psc;
        spm_write_vol(beta_img, beta_array);

    end
end