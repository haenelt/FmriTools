% VASO realignment plots

% This script computes the motion parameter plots for realigned vaso time
% series. BOLD and VASO motion parameters are shown in the same plot. This
% helps the identification of inconsistent bold and vaso motion
% corrections.

% created by Daniel Haenelt
% Date created: 20-10-2019
% Last modified: 25-05-2020

% array of of input time series
rp_bold = {
    '/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_1/data.nii',...
    };

rp_vaso = {
    '/data/pt_01880/Experiment1_ODC/p5/odc/VASO2/Run_1/data.nii',...
    };

% parameters
path_output = '';

%%% do not edit below %%%

% make folder
if ~exist(path_output,'dir') 
    mkdir(path_output); 
end

% plot motion regressors
for i = 1:length(rp_bold)
    
    % basename of output
    name_output = ['rp_ merge_run_' num2str(i)];
    
    % read realignment parameters
    M_bold = dlmread(rp_bold{i});
    M_vaso = dlmread(rp_vaso{i});
    
    % rad2deg
    M_bold(:,4) = radtodeg(M_bold(:,4));
    M_bold(:,5) = radtodeg(M_bold(:,5));
    M_bold(:,6) = radtodeg(M_bold(:,6));
    M_vaso(:,4) = radtodeg(M_vaso(:,4));
    M_vaso(:,5) = radtodeg(M_vaso(:,5));
    M_vaso(:,6) = radtodeg(M_vaso(:,6));
    
    % translational displacement
    M_bold_trans = sqrt(M_bold(:,1).^2+M_bold(:,2).^2+M_bold(:,3).^2);
    M_vaso_trans = sqrt(M_vaso(:,1).^2+M_vaso(:,2).^2+M_vaso(:,3).^2);
    
    % rotational displacement
    M_bold_rot = sqrt(M_bold(:,4).^2+M_bold(:,5).^2+M_bold(:,6).^2);
    M_vaso_rot = sqrt(M_vaso(:,4).^2+M_vaso(:,5).^2+M_vaso(:,6).^2);
    
    transFig = figure('visible','off');
    hold on
    plot(M_bold(:,1),'red');
    plot(M_bold(:,2),'blue');
    plot(M_bold(:,3),'green');
    plot(M_vaso(:,1),'red--');
    plot(M_vaso(:,2),'blue--');
    plot(M_vaso(:,3),'green--');
    plot(M_bold_trans,'black','LineWidth',2);
    plot(M_vaso_trans,'black--','LineWidth',2);
    title(['Translational movement in session ' name_output],'Interpreter','None');
    xlabel('number of volume');
    ylabel('Translation in mm');
    legend('x (BOLD)','y (BOLD)','z (BOLD)','x (vaso)','y (VASO)','z (VASO)','total translation (BOLD)','total translation (VASO)');
    saveas(gcf,fullfile(path_output, [name_output '_trans.png']));
    close(transFig);
    
    radFig = figure('visible','off');
    hold on
    plot(M_bold(:,4),'red');
    plot(M_bold(:,5),'blue');
    plot(M_bold(:,6),'green');
    plot(M_vaso(:,4),'red--');
    plot(M_vaso(:,5),'blue--');
    plot(M_vaso(:,6),'green--');
    plot(M_bold_rot,'black','LineWidth',2);
    plot(M_vaso_rot,'black--','LineWidth',2);
    title(['Rotational movement in session ' name_output],'Interpreter','None');
    xlabel('number of volume');
    ylabel('Rotation in deg');
    legend('pitch (BOLD)','roll (BOLD)','yaw (BOLD)','pitch (VASO)','roll (VASO)','yaw (VASO)','total rotation (BOLD)','total rotation (VASO)');
    saveas(gcf,fullfile(path_output,[name_output '_rot.png']));
    close(radFig);

end