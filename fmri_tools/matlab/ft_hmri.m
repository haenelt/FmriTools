function ft_hmri(file_t1w, file_pdw, file_b0, file_b1, pathHMRI_defaults, ...
    pathHMRI_b1_defaults)
% hMRI processing
%
% ft_hmri(file_t1w, file_pdw, file_b0, file_b1, pathHMRI_defaults, ...
%    pathHMRI_b1_defaults)
%
% Inputs:
%   file_t1w             - cell array of multi-echo t1w data set.
%   file_pdw             - cell array of multi-echo pdw data set.
%   file_b0              - cell array of B0 field map {magn1, magn2, phase diff.}
%   file_b1              - cell array of B1 field map.
%   pathHMRI_defaults    - filename of matlab file with hMRI default settings.
%   pathHMRI_b1_defaults - filename of matlab file with hMRI B1 default settings.
%   
% Computation of quantitative maps from acquired MPMs. The script is
% based on the batch GUI of the SPM12 toolbox. To run the script, both the
% hMRI and SPM12 have to be installed. The script gets the root directory 
% of the SPM12 toolbox from the matlab search path. All input files are 
% defined as cell arrays. The frame number in the filename string is 
% omitted since 1 (3d volume) is the default parameter for file selections. 
% The script consists of the following steps: (1) 7T configurations are 
% loaded. (2) Images are oriented to align with MNI space. Input for 
% registration is the first pd echo. The registration is then applied to 
% all other images (pdw, t1w, b1, b0). (3) Quantitative maps are created. 
% (4) Resulting maps are reoriented back to native space. (5) Maps are 
% scaled.
%
% Dependencies (should be installed as SPM toolbox):
% - hMRI toolbox: https://hmri-group.github.io/hMRI-toolbox/

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
