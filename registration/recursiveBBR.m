% recursive boundary-based registration

% This script performs recurive BBR to register a functional data set with
% surfaces generated from a separate whole-brain anatomy.

% parameter mode (explanation from one of Tim vanMourik's scripts:
% `r' for rotation around all axes, `rx', `ry', or `rz' for a single rotation around the respective 
% axis and for example `rxrz' for a rotation around the x-axis and the z-axis. In the same way the 
% scaling and translation can be modified. For example the combination `rystytz' would optimise for 
% a rotation around the y-axis, a scaling in all direction and a translation along the y-axis and 
% along the z-axis. The default is 'rst' for all transformations.

% created by Daniel Haenelt
% Date created: 09-12-2019
% Last modified: 09-12-2019

% input surfaces
input_white = {
    '/home/daniel/Schreibtisch/temp/lh.layer10_def_match',...
    '/home/daniel/Schreibtisch/temp/rh.layer10_def_match',...
    };

input_pial = {
    '/home/daniel/Schreibtisch/temp/lh.layer0_def_match',...
    '/home/daniel/Schreibtisch/temp/rh.layer0_def_match',...
    };

% parameters
subjects_dir = '/home/daniel/Schreibtisch';
ref_vol = 'temo/mean_data.nii';
mask = '';
min_vox = 50; % 10
min_vtx = 2000; % 500
cuboid = true;
tetrahedra = false;
nn_smooth = 0.5;
registration_mode = 'rst';
reverse_contrast = true; % if inner surface is darker than outer surface
cost_method = 'GreveFischl'; % GreveFischl, sum, sumOfSquares
contrast_method = 'gradient'; % gradient, average, fixedDistance, relativeDistance, extrema

% output filenames
output_surf = 'boundaries_out.mat';
output_cmap = 'sd.nii';

% path
pathSPM = '/home/daniel/source/spm12';
pathFSURF = '/usr/local/freesurfer';
pathMOURIK = '/home/daniel/projects/OpenFmriAnalysis';
pathLIB = '/home/daniel/projects/scripts/lib';

%%% do not edit below %%%

% file separator
fs = filesep;

% add paths to the interpreter's search path
addpath(genpath(pathSPM));
addpath(genpath(pathFSURF));
addpath(genpath(pathMOURIK));
addpath(genpath(pathLIB));

% parameters for bbr registration
cfg = [];
cfg.i_SubjectDirectory = subjects_dir;
cfg.i_ReferenceVolume = ref_vol;
cfg.i_Boundaries = 'boundaries_in.mat';
cfg.i_Mask = mask;
cfg.i_MinimumVoxels = min_vox;
cfg.i_MinimumVertices = min_vtx;
cfg.i_CuboidElements = cuboid;
cfg.i_Tetrahedra = tetrahedra;
cfg.i_NeighbourSmoothing = nn_smooth;
cfg.i_Mode = registration_mode;
cfg.i_ReverseContrast = reverse_contrast;
cfg.i_OptimisationMethod = cost_method;
cfg.i_ContrastMethod = contrast_method;
cfg.o_Boundaries = output_surf;
cfg.o_DisplacementMap = output_cmap;

% surf2mat
save_output = [cfg.i_SubjectDirectory fs cfg.i_Boundaries];
[wSurface, pSurface, faceData] = surf2mat(input_white, input_pial, pathFSURF, save_output);

% registration
tvm_recursiveBoundaryRegistration(cfg);

% load boundary output
load([cfg.i_SubjectDirectory fs cfg.o_Boundaries]);

% write surfaces
hemi = {'lh','rh'};
for i = 1:length(hemi)
    write_surf([cfg.i_SubjectDirectory fs hemi{i} '.white_def'], wSurface{i}(:,1:3), faceData{i} + 1);
    write_surf([cfg.i_SubjectDirectory fs hemi{i} '.pial_def'], pSurface{i}(:,1:3), faceData{i} + 1);
end