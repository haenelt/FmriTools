function ft_deface(file_in, file_other, write_mask)
% MRI defacing
%
% ft_deface(file_in, file_other, write_mask)
%
% Inputs:
%   file_in    - file name of input image.
%   file_other - file names of other images.
%   write_mask - save binary.
%
% This function defaces an MRI image using the SPM function spm_deface and
% writes the defaced image into the subfolder deface. The same defacing
% mask can be applied to a list of other images (e.g. useful if the data
% set consists of multiple echoes).

% constants
rtol = 1e-5; % relative tolerance
atol = 1e-8; % absolute tolerance

% default parameter
if ~exist('file_other', 'var')
    file_other = [];
end

if ~exist('write_mask', 'var')
    write_mask = false;
end

[path, basename, ext] = fileparts(file_in);
path_out = fullfile(path, 'deface');
file_tmp = fullfile(path, ['anon_' basename ext]);
file_out = fullfile(path_out, ['anon_' basename ext]);

% make output folder
if ~exist(path_out,'dir') 
    mkdir(path_out);
end

% deface input image
spm_deface(file_in);
movefile(file_tmp, file_out);

% get mask
img_in = spm_vol(file_in);
img_out = spm_vol(file_out);

arr_in = spm_read_vols(img_in);
arr_out = spm_read_vols(img_out);
arr_mask = abs(arr_in-arr_out) <= atol+rtol*abs(arr_out);

if write_mask
    img_out.fname = fullfile(path_out, ['mask' ext]);
    spm_write_vol(img_out, arr_mask);
end

% apply to all images
for i = 1:length(file_other)
    [~, basename, ext] = fileparts(file_other(i));
    file_out = fullfile(path_out, ['anon_' basename ext]);
    img_in = spm_vol(file_other(i));
    arr_in = spm_read_vols(img_in);
    arr_out = arr_in .* arr_mask;
    img_in.fname = file_out;
    spm_write_vol(img_in, arr_out);
end
