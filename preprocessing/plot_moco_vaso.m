% VASO realignment plots

% This script computes the motion parameter plots for realigned vaso time
% series. BOLD and VASO motion parameters are shown in the same plot. This
% helps the identification of inconsistent bold and vaso motion
% corrections.

% created by Daniel Haenelt
% Date created: 20-10-2019
% Last modified: 20-10-2019

% array of of input time series
img_input = {
    '/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_1/data.nii',...
    '/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_2/data.nii',...
    '/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_3/data.nii',...
    '/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_4/data.nii',...
    '/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_5/data.nii',...
    '/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_6/data.nii',...
    '/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_7/data.nii',...
    '/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_8/data.nii',...
    '/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_9/data.nii',...
    '/data/pt_01880/Experiment1_ODC/p5/odc/VASO1/Run_10/data.nii',...
    };

% basenames
basename_bold = 'bold';
basename_vaso = 'vaso';

%%% do not edit below %%%

% preprocessing summary
if length(img_input) > 1
    path_diagnosis = fullfile(fileparts(fileparts(img_input{1})),'diagnosis');
else
    path_diagnosis = fullfile(fileparts(img_input{1}),'diagnosis');
end

if ~exist(path_diagnosis,'dir') 
    mkdir(path_diagnosis); 
end

% plot motion regressors
for i = 1:length(img_input)
    
    % read realignment parameters
    M_bold = dlmread(fullfile(fileparts(img_input{i}),['rp_' basename_bold '.txt']));
    M_vaso = dlmread(fullfile(fileparts(img_input{i}),['rp_' basename_vaso '.txt']));
    
    transFig = figure('visible','off');
    hold on
    plot(M_bold(:,1),'red');
    plot(M_bold(:,2),'blue');
    plot(M_bold(:,3),'green');
    plot(M_vaso(:,1),'red--');
    plot(M_vaso(:,2),'blue--');
    plot(M_vaso(:,3),'green--');
    title(['Translational movement in session ' num2str(i)]);
    xlabel('number of volume');
    ylabel('Translation in mm');
    legend('x (bold)','y (bold)','z (bold)','x (vaso)','y (VASO)','z (VASO)');
    saveas(gcf,fullfile(path_diagnosis,['moco_mm_' basename_bold '_' basename_vaso '_' num2str(i) '.png']));
    close(transFig);
    
    radFig = figure('visible','off');
    hold on
    plot(M_bold(:,4),'red');
    plot(M_bold(:,5),'blue');
    plot(M_bold(:,6),'green');
    plot(M_vaso(:,4),'red--');
    plot(M_vaso(:,5),'blue--');
    plot(M_vaso(:,6),'green--');
    title(['Rotational movement in session ' num2str(i)]);
    xlabel('number of volume');
    ylabel('Rotation in rad');
    legend('pitch (BOLD)','roll (BOLD)','yaw (BOLD)','pitch (VASO)','roll (VASO)','yaw (VASO)');
    saveas(gcf,fullfile(path_diagnosis,['moco_rad_' basename_bold '_' basename_vaso '_' num2str(i) '.png']));
    close(radFig);

end