function ft_rBBR(input)
% Run rBBR
%
% ft_rBBR(input)
%
% Inputs:
%   input - mat-file with all input variables.
%
% This function calls the recursive BBR function in the OpenFmriAnalysis
% toolbox. From an loaded mat-file, all input parameters are ordered in a
% cell structure, i.e. all the needed variables should be contained in the
% input file.

% load input mat file
load(input);

% parameters for bbr registration
cfg = [];
cfg.i_SubjectDirectory = subjects_dir;
cfg.i_ReferenceVolume = reference_volume;
cfg.i_Boundaries = input_surf;
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

% registration
tvm_recursiveBoundaryRegistration(cfg);
