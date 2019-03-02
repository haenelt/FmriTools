function get_biopac_regressor(input, biopac_input, pathSPM, path_output, TR)
% This function computes nuisance regressors from peripheral cardiac and
% respiratory data (biopac).
% Inputs:
    % input: input time series
    % biopac_inpout: biopac *.mat file.
    % pathSPM: path to spm toolbox.
    % path_output: path where output is saved.
    % TR: repetition time in s.

% created by Daniel Haenelt
% Date created: 02-03-2019
% Last modified: 02-03-2019

% add paths to the interpreter's search path
addpath(pathSPM);
spm('defaults','FMRI');
spm_get_defaults('stats.maxmem',2^35); % maxmen indicates how much memory can be used
spm_get_defaults('cmdline',true); % no gui

% make output folder
if ~exist(path_output,'dir') 
    mkdir(path_output); 
end

% get paths and filenames
path_img = fileparts(input);
path_biopac = fileparts(biopac_input);

% convert biopac file to a textfile
biopac = load(biopac_input);
respiratory_data = biopac.data(:,1);
pulse_data = biopac.data(:,2);
save(fullfile(path_biopac,'respiratory_data.txt'),'respiratory_data','-ascii');
save(fullfile(path_biopac,'pulse_data.txt'),'pulse_data','-ascii');

% get time of first trigger
trigger_data = biopac.data(:,4);
exit = 0;
i = 0;
while ~exit
    i = i + 1;
    if trigger_data(i) ~= 0
        exit = 1;
    end
end

shift_log = i * biopac.isi / 1000;

% number of volumes
data_img = spm_vol(input);
nt = length(data_img);

% physio
matlabbatch{1}.spm.tools.physio.save_dir = {path_output};
matlabbatch{1}.spm.tools.physio.log_files.vendor = 'Custom';
matlabbatch{1}.spm.tools.physio.log_files.cardiac = {fullfile(path_biopac,'pulse_data.txt')};
matlabbatch{1}.spm.tools.physio.log_files.respiration = {fullfile(path_biopac,'respiratory_data.txt')};
matlabbatch{1}.spm.tools.physio.log_files.scan_timing = {''};
matlabbatch{1}.spm.tools.physio.log_files.sampling_interval = biopac.isi / 1000;
matlabbatch{1}.spm.tools.physio.log_files.relative_start_acquisition = shift_log;
matlabbatch{1}.spm.tools.physio.log_files.align_scan = 'first';
matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nslices = data_img(1).dim(3);
matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = [];
matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.TR = TR;
matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Ndummies = 0;
matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nscans = nt;
matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.onset_slice = round(data_img(1).dim(3)/2);
matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.time_slice_to_slice = [];
matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nprep = [];
matlabbatch{1}.spm.tools.physio.scan_timing.sync.nominal = struct([]);
matlabbatch{1}.spm.tools.physio.preproc.cardiac.modality = 'PPU';
matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = 'initial_cpulse_kRpeakfile.mat';
matlabbatch{1}.spm.tools.physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
matlabbatch{1}.spm.tools.physio.model.output_multiple_regressors = 'nuisance_regressor.txt';
matlabbatch{1}.spm.tools.physio.model.output_physio = 'physio.mat';
matlabbatch{1}.spm.tools.physio.model.orthogonalise = 'none';
matlabbatch{1}.spm.tools.physio.model.censor_unreliable_recording_intervals = false;
matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.c = 3;
matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.r = 4;
matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.cr = 1;
matlabbatch{1}.spm.tools.physio.model.rvt.no = struct([]);
matlabbatch{1}.spm.tools.physio.model.hrv.no = struct([]);
matlabbatch{1}.spm.tools.physio.model.noise_rois.no = struct([]);
matlabbatch{1}.spm.tools.physio.model.movement.no = struct([]);
matlabbatch{1}.spm.tools.physio.model.other.no = struct([]);
matlabbatch{1}.spm.tools.physio.verbose.level = 2;
matlabbatch{1}.spm.tools.physio.verbose.fig_output_file = '';
matlabbatch{1}.spm.tools.physio.verbose.use_tabs = false;

% run
spm_jobman('run',matlabbatch);

% clear matlabbatch
clear matlabbatch