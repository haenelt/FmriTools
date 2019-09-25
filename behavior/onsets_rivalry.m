% Compute rivalry onset data
% 
% From the condition logfile of the binocular rivalry paradigm, onsets and
% durations for the following GLM analysis are computed. Onsets are split
% into rest, left eye, right eye or piecemeal periods. Double button
% presses are defined as same button presses falling within a defined time
% window. Onsets are further excluded if they are not neighboured by
% correct switches (i.e., if they are neighboured by the same button
% presses which are above the double button press window) and if their
% durations (before and after onset) are below a set duration threshold.
% Correct switches also include that onsets are not neighboured by
% piecemeal periods.
%
% In the experiment, participants wear the red filter over their left eye
% and the green filter over their right eye. In SwitchData, button 
% responses of 2 and 1 refer to red and green onsets, respectively.
%
% created by Daniel Haenelt
% Date created: 05-08-2019         
% Last modified: 25-09-2019  

% input of condition files
input = {
    '/data/pt_01880/Experiment2_Rivalry/p4/rivalry/GE_EPI2/Run_1/logfiles/p4_GE_EPI2_Run1_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p4/rivalry/GE_EPI2/Run_2/logfiles/p4_GE_EPI2_Run2_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p4/rivalry/GE_EPI2/Run_3/logfiles/p4_GE_EPI2_Run3_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p4/rivalry/GE_EPI2/Run_4/logfiles/p4_GE_EPI2_Run4_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p4/rivalry/GE_EPI2/Run_5/logfiles/p4_GE_EPI2_Run5_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p4/rivalry/GE_EPI2/Run_6/logfiles/p4_GE_EPI2_Run6_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p4/rivalry/GE_EPI2/Run_7/logfiles/p4_GE_EPI2_Run7_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p4/rivalry/GE_EPI2/Run_8/logfiles/p4_GE_EPI2_Run8_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p4/rivalry/GE_EPI2/Run_9/logfiles/p4_GE_EPI2_Run9_rivalry_Cond.mat',...
    '/data/pt_01880/Experiment2_Rivalry/p4/rivalry/GE_EPI2/Run_10/logfiles/p4_GE_EPI2_Run10_rivalry_Cond.mat',...
    };

% parameters
response_threshold = 3; % in s

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
            % exclude if neighboured by piecemeal periods
            if data.SwitchData(j-1,1) == 3 || data.SwitchData(j+1,1) == 3
                data_consider(j) = 0;
            else
                data_consider(j) = 1;
            end
        else
            data_consider(j) = 0;
        end
    end

    % clean switches (only consider time points over time threshold)
    % outliers get the response 4 but are not written out
    response = data.SwitchData(:,1);
    response_time = data.SwitchData(:,2);
    response_duration = zeros(length(response),1);
    
    time_delta = response_time(2) - response_time(1);
    if time_delta < response_threshold
        response(1) = 4;
    end
    response_duration(1) = time_delta;
    
    for j = 2:length(response) - 1
        time_delta1 = response_time(j) - response_time(j-1);
        time_delta2 = response_time(j+1) - response_time(j);
        if time_delta1 < response_threshold || time_delta2 < response_threshold
            response(j) = 4;
        end
        response_duration(j) = time_delta2;
    end
    
    % clean responses
    response_time = response_time(data_consider == 1);
    response_duration = response_duration(data_consider == 1);
    response = response(data_consider == 1);
    
    response_time = response_time(response ~= 4);
    response_duration = response_duration(response ~= 4);
    response = response(response ~= 4);
    
    % save output as mat file
    names = cell(1,3);
    names{1} = 'rest';
    names{2} = 'left';
    names{3} = 'right';
    onsets = cell(1,3);
    onsets{1} = data.onsets{1}';
    onsets{2} = response_time(response==2);
    onsets{3} = response_time(response==1);
    durations = cell(1,3);
    durations{1} = [data.durations{1} ; data.durations{1}];
    durations{2} = response_duration(response==2);
    durations{3} = response_duration(response==1);
    
    % save finale condition file
    save(fullfile(pathfile,[name '_threshold' ext]),'names','onsets','durations');
    
    % plot design matrix
    nscale  = 10; % to shift the comma before rounding
    rest    = [ceil(nscale*data.onsets{1}') ceil(nscale*data.onsets{1}' + nscale*[data.durations{1} ; data.durations{1}])];
    left    = [ceil(nscale*response_time(response==2)) ceil(nscale*response_time(response==2) + nscale*response_duration(response==2))];
    right   = [ceil(nscale*response_time(response==1)) ceil(nscale*response_time(response==1) + nscale*response_duration(response==1))];
    mixed   = [ceil(nscale*response_time(response==3)) ceil(nscale*response_time(response==3) + nscale*response_duration(response==3))];
    
    M = zeros(1,1,3);
    for j = 1:length(rest(:,1))
        M(rest(j,1):rest(j,2),1,:) = 1;
    end
    if ~isempty(left)
        for j = 1:length(left(:,1))
            M(left(j,1):left(j,2),2,1) = 1;
        end
    end
    if ~isempty(right)
        for j = 1:length(right(:,1))
            M(right(j,1):right(j,2),3,2) = 1;
        end
    end
    if ~isempty(mixed)
        for j = 1:length(mixed(:,1))
            M(mixed(j,1):mixed(j,2),4,:) = 1;
        end
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