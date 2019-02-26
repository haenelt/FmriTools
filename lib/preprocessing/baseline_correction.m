function baseline_correction(input, TR, cutoff_highpass, pathSPM)

% This function computes a baseline correction of a functional time series.
% The output time series gets a prefix b to the file name.
% Inputs:
    % input: file name of time series.
    % TR: repetition time in s.
    % cutoff_highpass: highpass 1/cutoff frequency in Hz.
    % pathSPM: path to spm12 folder.

% created by Daniel Haenelt
% Date created: 07-12-2018
% Last modified: 10-12-2018

% add spm to path
addpath(pathSPM);

% get fileparts of input
[path, file, ext] = fileparts(input);

% load input time series
data_img = spm_vol(input);
data_array = spm_read_vols(data_img);

% get image dimensions
dim = data_img(1).dim;
nt = length(data_img);

% high-pass filtering of each voxel time-course
K.RT = TR;
K.row = 1:1:nt;
K.HParam = cutoff_highpass;

nK = spm_filter(K); % now we compute the low frequencies that should be removed from signal change
% nK.X0 contains the frequencies that can be removed. A second call of spm_filter removes these frequencies

data_array_corr = zeros([dim nt]);
for x = 1:dim(1)
    for y = 1:dim(2)
        for z = 1:dim(3)
            sig = data_array(x,y,z,:);
            sig = reshape(sig,[nt 1]);
            sig = spm_filter(nK,sig);
            sig = reshape(sig,[1 1 1 nt]);
            data_array_corr(x,y,z,:) = sig;
        end
    end
end

% write output
for i = 1:nt
  data_img(i).dim = dim;
  data_img(i).fname = fullfile(path, ['b' file ext]);
  spm_write_vol(data_img(i), data_array_corr(:,:,:,i));
end
