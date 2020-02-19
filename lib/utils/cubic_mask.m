function cubic_mask(input, path_output, name_output, c, r, pathSPM)
% This function computes a 3D cubic mask within and input array. The array
% dimensions are taken from the input nifti file. A binary mask is written.
% Inputs:
    % input: file name of time series.
    % path_output: path where output is written.
    % name_output: basename of output file.
    % c: array of center coordinates (x,y,z)
    % r: array of radius along coordinate axes (rx,ry,rz).
    % pathSPM: path to spm12 folder.

% created by Daniel Haenelt
% Date created: 19-02-2020
% Last modified: 19-02-2020

% add spm to path
addpath(pathSPM);

% make output folder
if ~exist(path_output,'dir') 
    mkdir(path_output); 
end

% get fileparts of input
[~, ~, ext] = fileparts(input);

% load input time series
data_img = spm_vol(input);

% get image dimensions
dim = data_img(1).dim;

% generate empty array
data_array = zeros(dim);

% get cube coordinates
c1 = c - r;
c2 = c + r;

% truncate cube at array borders
if c1(1) < 1
    c1(1) = 1;
end

if c1(2) < 1
    c1(2) = 1;
end

if c1(3) < 1
    c1(3) = 1;
end

if c2(1) > dim(1)
    c2(1) = dim(1);
end

if c2(2) > dim(2)
    c2(2) = dim(2);
end

if c2(3) > dim(3)
    c2(3) = dim(3);
end

% get binary cube
data_array(c1(1):c2(1),c1(2):c2(2),c1(3):c2(3)) = 1;

% write output
data_img(1).dim = dim;
data_img(1).fname = fullfile(path_output, [name_output ext]);
spm_write_vol(data_img(1), data_array);