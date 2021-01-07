function skullstrip_spm12(filename, path_output)
% Skullstrip SPM12.
%
% The computation of the skullstrip mask is done on the PD-weighted INV2 
% image. According to S. Kashyap, this shows the best results. Outputs are 
% written in a subfolder of the given output path.    

% created by Daniel Haenelt
% Date created: 01-11-2018             
% Last modified: 07-01-2021

% get SPM12 root directory from matlab search path
pathSPM = what('spm12').path;

% set spm default parameters
spm('defaults','FMRI');
spm_get_defaults('stats.maxmem',2^35); % maxmen indicates how much memory can be used
spm_get_defaults('cmdline',true); % no gui

% parameters
csf_max = 0.1; % c3 tissue class threshold
bone_max = 0.1; % c4 tissue class threshold
soft_max = 0.1; % c5 tissue class threshold
air_max = 0.1; % c6 tissue class threshold

% get path and basename
[~, basename, ext] = fileparts(filename);

% make output folder
path_skull = fullfile(path_output, 'skull');
if ~exist(path_skull, 'dir')
    mkdir(path_skull);
end

% copy input file to output folder
copyfile(filename, fullfile(path_skull, [basename ext]));

matlabbatch{1}.spm.spatial.preproc.channel.vols = {fullfile(path_skull, [basename ext])};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 18.0;
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {fullfile(pathSPM, 'tpm', 'TPM.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {fullfile(pathSPM, 'tpm', 'TPM.nii,2')};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {fullfile(pathSPM, 'tpm', 'TPM.nii,3')};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {fullfile(pathSPM, 'tpm', 'TPM.nii,4')};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {fullfile(pathSPM, 'tpm', 'TPM.nii,5')};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {fullfile(pathSPM, 'tpm', 'TPM.nii,6')};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0.0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.samp = 2.0;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];

% run
spm_jobman('run', matlabbatch);

% remove copied input file
delete(fullfile(path_skull, [basename ext]));

% load tissue classes
gm_img = spm_vol(fullfile(path_skull, ['c1' basename ext]));
wm_img = spm_vol(fullfile(path_skull, ['c2' basename ext]));
csf_img = spm_vol(fullfile(path_skull, ['c3' basename ext]));
bone_img = spm_vol(fullfile(path_skull, ['c4' basename ext]));
soft_img = spm_vol(fullfile(path_skull, ['c5' basename ext]));
air_img = spm_vol(fullfile(path_skull, ['c6' basename ext]));
    
gm_array = spm_read_vols(gm_img);
wm_array = spm_read_vols(wm_img);
csf_array = spm_read_vols(csf_img);
bone_array = spm_read_vols(bone_img);
soft_array = spm_read_vols(soft_img);
air_array = spm_read_vols(air_img);
    
% generate skullstrip mask
mask_array = gm_array + wm_array;
mask_array(mask_array > 0) = 1;
    
% get rid of pial noise
mask_array(csf_array >= csf_max) = 0;
mask_array(bone_array >= bone_max) = 0;
mask_array(soft_array >= soft_max) = 0;
mask_array(air_array >= air_max) = 0;
    
% save skullstrip mask
gm_img.fname = fullfile(path_skull, 'skullstrip_mask.nii');
spm_write_vol(gm_img, mask_array);
