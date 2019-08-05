% Compute rivalry onset data
% 
% From the condition logfile of the binocular rivalry paradigm, onsets and
% durations for the following GLM analysis are computed. Onsets are split
% into rest, left eye, right eye or piecemeal periods. Onsets are excluded
% if their durations are below a set time threshold. Optionally, an outlier
% regressor is written out.
%
% In the experiment, participant wear the red filter over their left eye
% and the green filter over their right eye. In SwitchData, button 
% responses of 2 and 1 refer to red and green onsets, respectively.
%
% created by Daniel Haenelt
% Date created: 05-08-2019         
% Last modified: 05-08-2019  

% input of condition files
input = {'/home/daniel/Schreibtisch/p1_GE_EPI1_Run1_rivalry_Cond.mat'};

% parameters
response_threshold = 4.5; % in s
write_outlier = true; % add outlier regressor

%%% do not edit below %%%

% loop through input files
for i = 1:length(input)

    % get filename
    [pathfile,name,ext] = fileparts(input{i});
    
    % load mat file
    data = load(input{i});
    
    % clean switches (delete repeating entries)
    data_repeat = zeros(length(data.SwitchData(:,1)),1);
    data_repeat(1) = 1;
    for j = 2:length(data.SwitchData(:,1))
        if data.SwitchData(j-1,1) == data.SwitchData(j,1)
            data_repeat(j) = 0;
        else
            data_repeat(j) = 1;
        end
    end
    
    % sum responses
    response = [4 ; data.SwitchData(data_repeat ~= 0,1)]; % 4 until first switch
    time1 = [data.onsets{2}(1) ; data.SwitchData(data_repeat ~= 0,2)];
    time2 = [time1(2:end) ; data.onsets{1}(2)];
    time_delta = time2 - time1;
    
    % threshold responses
    response(time_delta < response_threshold) = 4;
    
    % save output as mat file
    if write_outlier
        names = cell(1,5);
        names{1} = 'rest';
        names{2} = 'left';
        names{3} = 'right';
        names{4} = 'mixed';
        names{5} = 'outlier';
        onsets = cell(1,5);
        onsets{1} = data.onsets{1}';
        onsets{2} = time1(response==2);
        onsets{3} = time1(response==1);
        onsets{4} = time1(response==3);
        onsets{5} = time1(response==4);
        durations = cell(1,5);
        durations{1} = [data.durations{1} ; data.durations{1}];
        durations{2} = time_delta(response==2);
        durations{3} = time_delta(response==1);
        durations{4} = time_delta(response==3);
        durations{5} = time_delta(response==4);
    else
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
    end

    % save finale condition file
    save(fullfile(pathfile,[name '_threshold' ext]),'names','onsets','durations');
    
    % plot design matrix
    nscale  = 10; % to shift the comma before rounding
    rest    = [ceil(nscale*data.onsets{1}') ceil(nscale*data.onsets{1}' + nscale*[data.durations{1} ; data.durations{1}])];
    left    = [ceil(nscale*time1(response==2)) ceil(nscale*time1(response==2) + nscale*time_delta(response==2))];
    right   = [ceil(nscale*time1(response==1)) ceil(nscale*time1(response==1) + nscale*time_delta(response==1))];
    mixed   = [ceil(nscale*time1(response==3)) ceil(nscale*time1(response==3) + nscale*time_delta(response==3))];
    outlier = [ceil(nscale*time1(response==4)) ceil(nscale*time1(response==4) + nscale*time_delta(response==4))];
    
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
    for j = 1:length(outlier(:,1))
        M(outlier(j,1):outlier(j,2),5,:) = 1;
    end
    
    % plot design matrix
    figure(i)
    imagesc(M);
    names = {'baseline'; 'left'; 'right'; 'mixed'; 'outlier'};
    set(gca,'xtick',1:5,'xticklabel',names);
    yt = arrayfun(@num2str,get(gca,'ytick')/nscale,'un',0);
    set(gca,'yticklabel',yt);
    ylabel('time in TR');
    title(['Design matrix for run ' num2str(i)]);
    
end