function retino_fourier(input, freq, fix, period, TR, pathSPM)

% This function computes a fourier analysis for retinotopical data. Two
% time series are saved, namely the real and imaginary part of the Fourier
% transform. The data is thresholded by only including voxels above the 
% mean intensity. This is a very crude method for non-brain voxel
% exclusion. The output data gets the prefix r.
% Inputs:
    % input: file name of baseline corrected time series.
    % freq: number of cycles.
    % fix: pre and post run baseline block in s.
    % period: cycle period in s
    % TR: repetition time in s.
    % pathSPM: path to the SPM12 source folder.

% created by Daniel Haenelt
% Date created: 08-12-2018
% Last modified: 10-12-2018

% add spm to path
addpath(pathSPM);

% prepare path and file name
[path, file, ext] = fileparts(input);

% number of dropped volumes
drop_start = round(fix/TR) + round(mod(freq,1)*period/TR);
drop_end = round(fix/TR);
drop = sum([drop_start drop_end]);

% get time series header
data_img = spm_vol(input);

% number of volumes
nt = length(data_img);

% load time series
data_array = spm_read_vols(data_img);
data_array = data_array(:,:,:,drop_start+1:nt-drop_end);

% shift time to first dim in time series matrix
ts = shiftdim(data_array,3);

% normalise time series
vm = repmat(mean(ts,1), [nt - drop 1 1 1]); % voxel-wise mean
ts = (ts - vm) ./ vm * 100; % voxel-wise percent signal change

% fourier transform
ft = fft(ts); % Fourier transform

% phase shift to cf (-90 deg)
% we start the stimulus from the 3 o'clock position. However, the first
% quarter of the period is discarded. Therefore, we apply a phase shift to
% compensate for that
ft = ft * exp(-1i*mod(freq,1)*2*pi);

% Real and imaginary parts
cI = -imag(ft); % imaginary component (negative because we want cosine as 0)
cR = real(ft); % real component

% Mask the brain
mask = data_array(:,:,:,1);
threshold = mean(mask(:)); % intensity threshold
nbvox = find(mask < threshold); % non-brain voxels
% Remove non-brain voxels
for i = 1:size(cI,1)
    tempI = cI(i,:,:,:);
    tempI(nbvox) = NaN;
    tempR = cR(i,:,:,:);
    tempR(nbvox) = NaN;
    cI(i,:,:,:) = tempI;
    cR(i,:,:,:) = tempR;
end

% shift time to last dimension
cI = shiftdim(cI,1);
cR = shiftdim(cR,1);

% write output
ndata_img = data_img(1:size(cI,4));
for i = 1:size(cI,4)
    ndata_img(i).fname = fullfile(path, ['r' file '_imag' ext]);
    ndata_img(i).pinfo = [1 0 0]';
    ndata_img(i).dt(1) = 16; % convert data type to double (float32)
    spm_write_vol(ndata_img(i),cI(:,:,:,i));
end
disp(['Saved imaginary image: ' ndata_img(1).fname]);

ndata_img = data_img(1:size(cR,4));
for i = 1:size(cR,4)
    ndata_img(i).fname = fullfile(path, ['r' file '_real' ext]);
    ndata_img(i).pinfo = [1 0 0]';
    ndata_img(i).dt(1) = 16; % convert data type to double (float32)
    spm_write_vol(ndata_img(i),cR(:,:,:,i));
end
disp(['Saved real image: ' ndata_img(1).fname]);


