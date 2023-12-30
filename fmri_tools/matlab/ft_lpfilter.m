function ft_lpfilter(input, TR, cutoff_lowpass, order)
% Lowpass filter
%
% ft_lpfilter(input, TR, cutoff_lowpass, order)
%
% Inputs:
%   input          - file name of time series.
%   TR             - repetition time in s.
%   cutoff_lowpass - lowpass 1/cutoff frequency in Hz.
%   order          - order of butterworth filter.
%
% This function computes a butterworth lowpass filter and applies it to a
% functional time series. The output gets a prefix l to the file name.

% get fileparts of input
[path, file, ext] = fileparts(input);

% load input time series
data_img = spm_vol(input);
data_array = spm_read_vols(data_img);

% get image dimensions
dim = data_img(1).dim;
nt = length(data_img);

% filter parameters
f.fs = 1/TR; % sampling rate
f.fc = 1/cutoff_lowpass; % cutoff frequency
f.order = order; % filter order

% butterworth filter
[b,a] = butter(f.order,f.fc/(f.fs/2));

% filter time series
data_array_corr = zeros([dim nt]);
for x = 1:dim(1)
    for y = 1:dim(2)
        for z = 1:dim(3)
            sig = data_array(x,y,z,:);
            data_array_corr(x,y,z,:) = filter(b,a,sig);
        end
    end
end

% write output
for i = 1:nt
  data_img(i).dim = dim;
  data_img(i).fname = fullfile(path, ['l' file ext]);
  spm_write_vol(data_img(i), data_array_corr(:,:,:,i));
end
