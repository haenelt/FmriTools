cond_input = {
    '/data/pt_01880/Experiment6_Albinism/p1/odc/GE_EPI1/Run_1/logfiles/p1_GE_EPI1_Run1_odc_Cond.mat',...
    '/data/pt_01880/Experiment6_Albinism/p1/odc/GE_EPI1/Run_2/logfiles/p1_GE_EPI1_Run2_odc_Cond.mat',...
    '/data/pt_01880/Experiment6_Albinism/p1/odc/GE_EPI1/Run_3/logfiles/p1_GE_EPI1_Run3_odc_Cond.mat',...
    };

cond = load(cond_input{1});

nruns = 4;
nconds = 3;
hrf_derivative = true;

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

% get contrast vector
if hrf_derivative == true
    if nconds == 2
        c = [-1 0 1 0 ; 
              1 0 -1 0 ; 
              1 0 0 0 ; 
              0 0 1 0 ; 
              0 -1 0 1 ; 
              0 1 0 -1 ; 
              0 1 0 0 ; 
              0 0 0 1];
    elseif nconds == 3
        c = [0 0 -1 0 1 0 ; 0 0 1 0 -1 0 ; -1 0 0 0 1 0 ; -1 0 1 0 0 0 ; 1 0 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 0 1 0 ;
            0 0 0 -1 0 1 ; 0 0 0 1 0 -1 ; 0 -1 0 0 0 1 ; 0 -1 0 1 0 0 ; 0 1 0 0 0 0 ; 0 0 0 1 0 0 ; 0 0 0 0 0 1];
    elseif nconds == 4
        c = [1 0 -1 0 0 0 0 0 ; -1 0 1 0 0 0 0 0 ; 1 0 0 0 -1 0 0 0 ; 0 0 1 0 0 0 -1 0 ; 1 0 0 0 0 0 0 0 ; 0 0 1 0 0 0 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 1 0 ; 
            0 1 0 -1 0 0 0 0 ; 0 -1 0 1 0 0 0 0 ; 0 1 0 0 0 -1 0 0 ; 0 0 0 1 0 0 0 -1 ; 0 1 0 0 0 0 0 0 ; 0 0 0 1 0 0 0 0 ; 0 0 0 0 0 1 0 0 ; 0 0 0 0 0 0 0 1];
    end
else
    if nconds == 2
        c = [-1 1 ; 1 -1 ; 1 0 ; 0 1];
    elseif nconds == 3
        c = [0 -1 1 ; 0 1 -1 ; -1 0 1 ; -1 1 0 ; 1 0 0 ; 0 1 0 ; 0 0 1];
    elseif nconds == 4
        c = [1 -1 0 0 ; -1 1 0 0 ; 1 0 -1 0 ; 0 1 0 -1 ; 1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1];
    end
end


contrast = [];
for i=1:nruns
    contrast = [contrast c/nruns];
    %contrast5 = [contrast5 zeros(1,nconds)];
end

% get constant term
contrast = [contrast zeros(size(contrast,1),nruns)];
contrast = [contrast ; zeros(1,size(contrast,2))];
contrast(end,end+1-nruns:end) = 1/nruns;


% include also if for hrf derivative

matlabbatch{1}.spm.stats.con.spmmat = {fullfile('bla','SPM.mat')};
matlabbatch{1}.spm.stats.con.delete = 0; % append contrast to old contrasts

for i = 1:length(name_contrast)
    matlabbatch{1}.spm.stats.con.consess{i}.tcon.name = name_contrast{i};
    matlabbatch{1}.spm.stats.con.consess{i}.tcon.weights = contrast(i,:);
    matlabbatch{1}.spm.stats.con.consess{i}.tcon.sessrep = 'none';
end


%%

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
