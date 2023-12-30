% TDM analysis
%
% This script executes the TDM method developed by Kendrick Kay et al.
% (2019) to separate early and late BOLD contributions within an fmri
% timeseries. In this script, a simple block design with one experimental
% and one baseline condition is expected. Furthermore, the onset times and
% stimulus durations of the experimental conditions are assumed to be the
% same across functional runs.
%
% Dependencies (should be added to the matlab search path):
% - knkutils: https://github.com/kendrickkay/knkutils
% - GLMdenoise: https://github.com/kendrickkay/GLMdenoise
% - TDM: https://github.com/kendrickkay/TDM

% input data
img_input = {
    '/data/pt_01983/func/flicker/GE_EPI1/Run_1/uadata.nii',...
    '/data/pt_01983/func/flicker/GE_EPI1/Run_2/uadata.nii',...
    '/data/pt_01983/func/flicker/GE_EPI1/Run_3/uadata.nii',...
    '/data/pt_01983/func/flicker/GE_EPI1/Run_4/uadata.nii',...
    '/data/pt_01983/func/flicker/GE_EPI1/Run_5/uadata.nii',...
    };

cond_input = {
    '/data/pt_01983/func/flicker/GE_EPI1/Run_1/logfiles/p1_GE_EPI1_Run1_flicker_Cond.mat',...
    '/data/pt_01983/func/flicker/GE_EPI1/Run_2/logfiles/p1_GE_EPI1_Run2_flicker_Cond.mat',...
    '/data/pt_01983/func/flicker/GE_EPI1/Run_3/logfiles/p1_GE_EPI1_Run3_flicker_Cond.mat',...
    '/data/pt_01983/func/flicker/GE_EPI1/Run_4/logfiles/p1_GE_EPI1_Run4_flicker_Cond.mat',...
    '/data/pt_01983/func/flicker/GE_EPI1/Run_5/logfiles/p1_GE_EPI1_Run5_flicker_Cond.mat',...
    };

% parameters
TR = 1; % repetition time in s
stim_duration = 30; % stimulus duration in s
sigma = 2; % sigma for gaussian blurring (bias corrected epi)
mask_threshold = 0.8; % masking threshold for implicit mask
fir_steps = 60; % number of time steps for fir model
basename = 'on'; % basename of output files
output_folder = 'contrast_tdm'; % name of folder where glm output is saved

%%% do not edit below %%%

% output folder is taken from the first entry of the input list
if length(img_input) > 1
    path_output = fullfile(fileparts(fileparts(img_input{1})),output_folder);
    path_psc = fullfile(fileparts(fileparts(path_output)),'results','psc','native');
else
    path_output = fullfile(fileparts(img_input{1}),output_folder);
    path_psc = fullfile(fileparts(path_output),'results','psc','native');
end

if ~exist(path_output,'dir')
    mkdir(path_output);
end

if ~exist(path_psc,'dir')
    mkdir(path_psc);
end

% get data
data = get_TDMdata(img_input);

% get design
% We take the condition file of the first run assuming that all runs have
% the same condition onsets and duratations. The condition file contains
% two conditions (flicker, baseline) and we only extract the onsets of the
% flicker condition.
design = get_TDMdesign(cond_input{1}, TR, size(data{1},2), size(data,2));

% get epi
% The corrected epi intensities are computed from the first run.
bc = get_TDMepi(img_input{1}, sigma);

% get mask
mask = get_TDMmask(data, mask_threshold);

% mask data
bc = bc(mask==1);
for i = 1:length(data)
    data{i} = data{i}(mask==1,:);
end

% get lookup table for reshape
imglookup = linspace(1,length(mask),length(mask));
imglookup = imglookup(mask==1);

% fit FIR model
% Dervive unconstrained estimates of response timecourses by fitting a
% finite impulse response (FIR) model to the data. fir_steps is the number
% of time points (trial onset plus fir_steps); these time points come at a
% rate of TR.
resultsFIR = GLMdenoisedata(...
    design,...
    data,...
    stim_duration,...
    TR,...
    'fir',...
    fir_steps,...
    struct('numboots',0,'numpcstotry',0),... % no bootstrapping, no denoising
    fullfile(path_output,'GLMdenoise_FIR_figures')...
    );

% determine an R2 threshold for the FIR results
% We want to collect timecourses from voxels that have some minimum
% level of signal-to-noise ratio for experimentally evoked BOLD responses.
% Use an automatic method for determining a reasonable threshold.
r2thresh = findtailthreshold(resultsFIR.R2(:));

% extract FIR timecourses
ix = resultsFIR.R2 > r2thresh; % voxels that pass the threshold
timecourses = resultsFIR.modelmd(ix,:,:); % N x fir_step + 1

% Get bias-corrected EPI intensities for the same voxels
bcvalues = bc(ix);

% perform TDM
% Now we will pass the timecourses to the core of the TDM technique, i.e.,
% to derive estimate of the latent early and late timecourses. In the TDM
% output figures, black and gray dots indicate the early and late response,
% respectively.
close all % close all open figures
resultsTDM = extracthrfmanifold(...
    permute(timecourses, [1 3 2]),...
    bcvalues,...
    TR,...
    {fullfile(path_output,'TDM_figures') -1},... % -1 indicates to also write eps figures
    struct()... % use defaults for everything
    );

% fit TDM GLM
% Now that we have derived Early and Late timecourses, we can use them to
% re-fit the fMRI time-series data.

% (1) perform convolution with the early and late timecourses separately
design_conv = design;
for p = 1:length(design_conv)
  temp = [];
  for q = 1:2
      temp = [temp conv2(full(design_conv{p}),resultsTDM.elhrf(q,:)')];
  end
  design_conv{p} = temp(1:size(design_conv{p},1),:);
end

% plot the design for the first run for sanity check. The first and the
% second column corresponds to the early and late timecourse, respectively.
figure;
imagesc(design_conv{1});
colormap(gray);

% fit GLM with TDM-derived timecourses
% Note that we already convolved in the HRFs, so we use an <hrfknobs> input of 1 which
% effectively convolves our design matrices with 1 (which doesn't change anything).
resultsTDMGLM = GLMdenoisedata(...
    design_conv,...
    data,...
    stim_duration,...
    TR,...
    'assume',...
    1,...
    struct('numboots',0,'numpcstotry',0),... % no bootstrapping, no denoising
    fullfile(path_output,'GLMdenoise_TDMGLM_figures')...
    );

% betas for early and late response of the on condition
betas0 = resultsTDMGLM.modelmd{2};

% reshape and write output
data_header = spm_vol(img_input{1});
data_header = data_header(1);

nx = data_header.dim(1);
ny = data_header.dim(2);
nz = data_header.dim(3);
nvox = nx * ny * nz;

name_file = {'tdm_early', 'tdm_late'};
for i = 1:length(name_file)
        
    % get new filename
    data_header.fname = fullfile(path_psc, ['psc_' basename '_' name_file{i} '.nii']);

    % reshape output
    betas0_out = zeros(nvox,1);
    betas0_out(imglookup) = betas0(:,i);
    betas0_out = reshape(betas0_out, nx, ny, nz);

    % write output
    spm_write_vol(data_header, betas0_out);
    
end



function [data] = get_TDMdata(file_in)

% This function sorts timeseries data into a cell array. Each timeseries is
% reshaped into a 2D array (ind x time).
% Inputs:
    % file_in: cell array of filenames.
% Outputs:
    % data: cell array of data in kendrick kay's format.

data = {};
for i = 1:length(file_in)
    
    % load data
    temp = spm_vol(file_in{i});
    temp_array = spm_read_vols(temp);
    temp_array(isnan(temp_array)) = 0;
    
    % reshape each time series
    data{i} = reshape(temp_array, [], length(temp));
end

end


function [design] = get_TDMdesign(mat_in, TR, nvols, nruns)

% This function converts the experimental conditions from the condition
% file in spm12 format to a design matrix in kendrick kay's convention
% for use in TDM. We only use an on-off block design and only convert the
% onsets of the on blocks to the design matrix.
% Inputs:
    % mat_in: mat-file containing condition information in spm12 format.
    % TR: TR in s.
    % nvols: number of volumes.
    % nruns: number of within-session runs.
% Outputs:
    % design: design matrix in kendrick kay's format.

% load condition file
design_mat = load(mat_in);

% get design matrix for tdm
nconds = 1; % only consider on blocks which is the first condition in our case!
design = {};
for i = 1:nruns
    design{i} = zeros(TR*nvols, nconds);
    for j = 1:nconds
        design{i}(round(design_mat.onsets{j}),j) = 1;
    end
    design{i} = sparse(design{i});
end

end


function [epi_reshape] = get_TDMepi(file_in, sigma, path_output, write_output)

% This function computes relative epi intensities by taking the ratio of
% the temporal mean and the gaussian filtered temporal mean. 
% Inputs:
    % file_in: filename of input timeseries.
    % sigma: sigma for gaussian filtering.
    % path_output: path where output is written.
    % write_output: write output files.
% Outputs:
    % epi_reshape: relative epi intensities sorted into 1D vector.

if nargin < 4
    path_output = '';
    write_output = false;
end

if nargin < 3
    sigma = 2;
end

% get filename
[~, name_file, ext_file] = fileparts(file_in);

% load data 
data = spm_vol(file_in);

% set nan to zero
data_array = spm_read_vols(data);
data_array(isnan(data_array)) = 0;

% get time average
data_array_mean = mean(data_array, 4);

% filter time average by gaussian
data_array_gaussian = imgaussfilt3(data_array_mean, sigma);

% get ratio
data_array_res = data_array_mean ./ data_array_gaussian;

% reshape to 1D vector
epi_reshape = reshape(data_array_res, [], 1);

% write output
if write_output == true

    data = data(1);
    data.fname = fullfile(path_output, [name_file '_vein' ext_file]);
    spm_write_vol(data, data_array_res);

    data.fname = fullfile(path_output, [name_file '_gauss' ext_file]);
    spm_write_vol(data, data_array_gaussian);

end

end


function [mask_array] = get_TDMmask(data, mask_threshold)

% This function calculates an implicit mask for GLM.
% Inputs:
    % data: data cell array in kendrick kay's format.
    % mask_threshold: threshold for binary mask.
% Outputs:
    % mask_array: binary mask of data vector.

if nargin < 2
    mask_threshold = 0.8;
end

% get dat mean
nruns = length(data);
data_mean = 0;
mask_array = zeros(size(data{1},1),1);
for i = 1:nruns
    data_mean = data_mean + mean(data{i}(:));
    mask_array = mask_array + mean(data{i},2);
end
data_mean = data_mean / nruns;
mask_array = mask_array / nruns;

% get binary mask
mask_array(mask_array < mask_threshold*data_mean) = 0;
mask_array(mask_array ~= 0) = 1;

end
