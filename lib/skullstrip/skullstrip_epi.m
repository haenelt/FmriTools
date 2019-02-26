function mask_array = skullstrip_epi(input, pathSPM, roi_size, scale, nerode, ndilate, savemask)
% Skullstrip input volume by defining an intensity threshold from the inner 
% of the brain volume. From a defined mid-point, a brain mask is grown 
% inside the brain. A binary filling holes algorithm is applied. To reduce 
% remaining skull within the brain mask, the mask is eroded and dilated 
% several times.    
% Inputs:
%   *input: input file.
%   *pathSPM: path to spm toolbox.
%   *roi_size: size of cubic roi for image intensity threshold
%   *scale: scale image intensity threshold
%   *nerode: number of eroding iterations
%   *ndilate: number of dilating iterations
%   *savemask: save mask time series (boolean).
%   *cleanup: delete intermediate files after running.
%    
% created by Daniel Haenelt
% Date created: 06-11-2018             
% Last modified: 25-02-2019

% default parameters
if ~exist('roi_size','var')
    roi_size = 10;
end

if ~exist('scale','var')
    scale = 0.7;
end

if ~exist('nerode','var')
    nerode = 3;
end

if ~exist('ndilate','var')
    ndilate = 1;
end

if ~exist('savemask','var')
    savemask = true;
end

% add spm toolbox
addpath(pathSPM);

% prepare path and filename
[path, file, ext] = fileparts(input);

% load data
data_img = spm_vol(input);
data_array = spm_read_vols(data_img);

% calculate mean intensity
data_mean = mean(data_array(:));

% get point within the brain
[x, y, z] = ind2sub(size(data_array),find(data_array > data_mean));
x_mean = round((max(x)+min(x))/2);
y_mean = round((max(y)+min(y))/2);
z_mean = round((max(z)+min(z))/2);

% initialize mask
mask_array = zeros(size(data_array));
mask_temp_array = zeros(size(data_array));
mask_array(x_mean,y_mean,z_mean) = 1;
mask_temp_array(x_mean,y_mean,z_mean) = 1;

% compute threshold
roi = data_array(...
    ceil(x_mean-roi_size/2):floor(x_mean+(roi_size-1)/2),...
    ceil(y_mean-roi_size/2):floor(y_mean+(roi_size-1)/2),...
    ceil(z_mean-roi_size/2):floor(z_mean+(roi_size-1)/2)...
    );
roi_mean = mean(roi(:));

% grow mask
while true
    
    % get ROI coordinates
    [x, y, z] = ind2sub(size(mask_array),find(mask_array > 0));
    coords = [x y z];
    
    % edge
    coords = coords(coords(:,1) > 1,:);
    coords = coords(coords(:,2) > 1,:);
    coords = coords(coords(:,3) > 1,:);
    coords = coords(coords(:,1) < size(data_array,1)-1,:);
    coords = coords(coords(:,2) < size(data_array,2)-1,:);
    coords = coords(coords(:,3) < size(data_array,3)-1,:);
        
    % calculate nearest neighbors
    coords1 = [coords(:,1)-1, coords(:,2), coords(:,3)];
    coords2 = [coords(:,1), coords(:,2)-1, coords(:,3)];
    coords3 = [coords(:,1), coords(:,2), coords(:,3)-1];
    coords4 = [coords(:,1)+1, coords(:,2), coords(:,3)];
    coords5 = [coords(:,1), coords(:,2)+1, coords(:,3)];
    coords6 = [coords(:,1), coords(:,2), coords(:,3)+1];
    
    % merge new coordinates
    coords = [coords ; coords1 ; coords2 ; coords3 ; coords4 ; coords5 ; coords6];
    
    % convert to indices
    inds = sub2ind(size(mask_array),coords(:,1),coords(:,2),coords(:,3));
    
    % get temporary mask from coordinates
    mask_temp_array(inds) = 1;
       
    % delete all old mask elements
    mask_temp_array(mask_array == 1) = 0;
    mask_temp_array(data_array < scale*roi_mean) = 0;
       
    % reinitialize mask temp
    mask_array = mask_array + mask_temp_array;
       
    % check break condition
    if length(find(mask_temp_array > 0)) == 0
        break
    end
    
end

% flood filling brain mask
mask_array = imfill(mask_array,'holes');

% erode mask
for i = 1:nerode
    mask_array = imerode(mask_array,ones(3));
end

% dilate mask
for i = 1:ndilate
    mask_array = imdilate(mask_array, ones(3));
end
    
% write maskedimage
data_array = data_array .* mask_array;
data_img.fname = fullfile(path, ['p' file ext]);
spm_write_vol(data_img, data_array);

% write mask
if savemask
    data_img.fname = fullfile(path, ['mask' file ext]);
    spm_write_vol(data_img, mask_array);
end