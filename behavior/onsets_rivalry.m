% Compute rivalry onset data
% 
% From the condition logfile of the binocular rivalry paradigm, onsets and
% durations for the following GLM analysis are computed. Onsets are split
% into rest, left eye, right eye or piecemeal periods. Double button
% presses are defined as same button presses falling within a defined time
% window. Onsets are further excluded if their durations are below a set
% time threshold and if they are not neighboured by correct switches (i.e.,
% if they are neighboured by the same button presses which are above the
% double button press window). Optionally, a mixed (save_mixed) regressor
% is written out.
%
% In the experiment, participants wear the red filter over their left eye
% and the green filter over their right eye. In SwitchData, button 
% responses of 2 and 1 refer to red and green onsets, respectively.
%
% created by Daniel Haenelt
% Date created: 05-08-2019         
% Last modified: 15-09-2019  

% input of condition files
input = {
    '/data/pt_01880/Experiment2_Rivalry/p1/rivalry/GE_EPI2/Run_1/logfiles/p1_GE_EPI2_Run1_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p1/rivalry/GE_EPI2/Run_2/logfiles/p1_GE_EPI2_Run2_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p1/rivalry/GE_EPI2/Run_3/logfiles/p1_GE_EPI2_Run3_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p1/rivalry/GE_EPI2/Run_4/logfiles/p1_GE_EPI2_Run4_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p1/rivalry/GE_EPI2/Run_5/logfiles/p1_GE_EPI2_Run5_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p1/rivalry/GE_EPI2/Run_6/logfiles/p1_GE_EPI2_Run6_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p1/rivalry/GE_EPI2/Run_7/logfiles/p1_GE_EPI2_Run7_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p1/rivalry/GE_EPI2/Run_8/logfiles/p1_GE_EPI2_Run8_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p1/rivalry/GE_EPI2/Run_9/logfiles/p1_GE_EPI2_Run9_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p1/rivalry/GE_EPI2/Run_10/logfiles/p1_GE_EPI2_Run10_rivalry_Cond.mat',...
    };

% parameters
response_threshold = 3; % in s
save_mixed = false;

%%% do not edit below %%%

% loop through input files
for i = 1:length(input)

    % get filename
    [pathfile,name,ext] = fileparts(input{i});
    
    % load mat file
    data = load(input{i});
    
    % remove double button presses (define double button press if the same
    % button is pressed within a defined time window)
    double_press_threshold = 0.2;
    for j = 1:length(data.SwitchData(:,1)) - 1
        if data.SwitchData(j,1) == data.SwitchData(j+1,1)
            time_delta = data.SwitchData(j+1,2) - data.SwitchData(j,2);
            if time_delta < double_press_threshold
                data.SwitchData(j+1,:) = [];
            end
        end
    end
    
    % clean switches (only consider time points with correct neighbouring)
    data_consider = zeros(length(data.SwitchData(:,1)),1);
    if data.SwitchData(1,1) ~= data.SwitchData(2,1) 
        data_consider(1) = 1;
    end
    for j = 2:length(data.SwitchData(:,1)) - 1
        if data.SwitchData(j-1,1) ~= data.SwitchData(j,1) && data.SwitchData(j,1) ~= data.SwitchData(j+1,1)
            data_consider(j) = 1;
        else
            data_consider(j) = 0;
        end
    end
    
    % sum responses
    response = data.SwitchData(data_consider ~= 0,1);
    time1 = data.SwitchData(data_consider ~= 0,2);
    time2 = data.SwitchData(find(data_consider ~= 0) + 1,2);
    time_delta = time2 - time1;
    
    % threshold responses (outliers get the response 4 but are not written
    % out)
    response(time_delta < response_threshold) = 4;
    
    % save output as mat file
    if save_mixed
        names = cell(1,4);
        names{1} = 'rest';
        names{2} = 'left';
        names{3} = 'right';
        names{4} = 'mixed';
        onsets = cell(1,4);
        onsets{1} = data.onsets{1}';
        onsets{2} = time1(response==2);
        onsets{3} = time1(response==1);
        onsets{4} = time1(response==3);
        durations = cell(1,4);
        durations{1} = [data.durations{1} ; data.durations{1}];
        durations{2} = time_delta(response==2);
        durations{3} = time_delta(response==1);
        durations{4} = time_delta(response==3);
    else
        names = cell(1,3);
        names{1} = 'rest';
        names{2} = 'left';
        names{3} = 'right';
        onsets = cell(1,3);
        onsets{1} = data.onsets{1}';
        onsets{2} = time1(response==2);
        onsets{3} = time1(response==1);
        durations = cell(1,3);
        durations{1} = [data.durations{1} ; data.durations{1}];
        durations{2} = time_delta(response==2);
        durations{3} = time_delta(response==1);
    end
    
    % save finale condition file
    save(fullfile(pathfile,[name '_threshold' ext]),'names','onsets','durations');
    
    % plot design matrix
    nscale  = 10; % to shift the comma before rounding
    rest    = [ceil(nscale*data.onsets{1}') ceil(nscale*data.onsets{1}' + nscale*[data.durations{1} ; data.durations{1}])];
    left    = [ceil(nscale*time1(response==2)) ceil(nscale*time1(response==2) + nscale*time_delta(response==2))];
    right   = [ceil(nscale*time1(response==1)) ceil(nscale*time1(response==1) + nscale*time_delta(response==1))];
    mixed   = [ceil(nscale*time1(response==3)) ceil(nscale*time1(response==3) + nscale*time_delta(response==3))];
    
    M = zeros(1,1,3);
    for j = 1:length(rest(:,1))
        M(rest(j,1):rest(j,2),1,:) = 1;
    end
    for j = 1:length(left(:,1))
        M(left(j,1):left(j,2),2,1) = 1;
    end
    for j = 1:length(right(:,1))
        M(right(j,1):right(j,2),3,2) = 1;
    end
    for j = 1:length(mixed(:,1))
        M(mixed(j,1):mixed(j,2),4,:) = 1;
    end
    
    % plot design matrix
    figure(i)
    set(gca,'visible','off');
    imagesc(M);
    names = {'baseline'; 'left'; 'right'; 'mixed'};
    set(gca,'xtick',1:4,'xticklabel',names);
    yt = arrayfun(@num2str,get(gca,'ytick')/nscale,'un',0);
    set(gca,'yticklabel',yt);
    ylabel('time in TR');
    title(['Design matrix for run ' num2str(i)]);
    saveas(gca,fullfile(pathfile,[name '_plot.png']));
    
end