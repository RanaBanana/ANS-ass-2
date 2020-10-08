% Load assignment2_data.mat

% Combine the events_ts and events_type to mark when the stimulus was
% presented.
table_ts_type = table(events_ts, events_type);

% Set of variables. Set to microseconds to match data.
timeBin = 25000; % time bin in microseconds
preStim = 500000; % time before stimulus onset in microseconds
postStim = 1000000; % time after stimulus offset in microseconds

% Extract values from table corresponding to the on- and offset
onTimes = table_ts_type.events_ts(events_type==1);
offTimes = table_ts_type.events_ts(events_type==31);

% Compute timeframe of PTSH (-500stim+1000ms). timeframePSTH is a table
% with two columns.
timeframePSTH = [onTimes-preStim, offTimes+postStim];
% spikes_PSTH = zeros(1, 19); % preallocation of the spikes_PSTH wouldn't work for 
% different sizes of arrays.

for i = 1:10 % actually for length(timeframePSTH), but your computer will crash and matlab will never forget. Matlab will haunt you till the end of days. So don't do this or you'll regret it for the rest of your life.
        spikes_PSTH = spikes_ts(1,[spikes_ts<timeframePSTH(i,2)]);
        % Make a histogram of the spikes within the PSTH (in each 25ms bin). 
        % It looks at the 500ms prestimulus to the 1000ms poststimulus 
        % interval, in steps of 'timeBin', which is 25 ms. 
        figure(i);
        histogram(spikes_PSTH,[timeframePSTH(i,1):timeBin:timeframePSTH(i,2)]);
        % Spice up the graph ~~*pretty pretty*~~ glam, ah, such wow
        grid on;
        xlabel('Time (microseconds)', 'FontSize', 12);
        ylabel('Spike Count', 'FontSize', 12);
        title('Histogram of Spikes per 25 ms time bin', 'FontSize', 14);
        
        % We probably have to use 'hold on' and 'hold off' to add different
        % timeframePTSHs in one graph.
end




%% Excercise 4
% Compute and plot the spike-triggered LFP signal (i.e. the average LFP [value of a] time signal
% occurring around the occurrence of an action potential) in the [-250, 250] ms window
% around the spike times. Do this separately for spikes in the baseline window and for
% spikes in the stimulus presentation window. What do you notice? Is there any difference
% in the LFP between baseline and stimulus presentation?

% Variable definitions
window_spike = 250000;
% lfp_data is expressed in uV.
LFP_table = table(lfp_data, lfp_ts);
%Test on first 12000 columns of LFP:
%small_data = LFP(:,1:12000);
%plot(small_data(2,:),small_data(1,:));

% Plotting the whole LFP signal; WARNING: ONLY RUN ON FAST COMPUTER
% plot_LFP = plot(LFP(2,:), LFP(1,:))



% New:

% Compute the 250 ms time window around the spike times (in columns).
window250_spikes = [spikes_ts-window_spike; spikes_ts+window_spike]';

% Preallocate the window250_lfp and the timewindow_lfp (both have two
% columns)
window250_lfp = zeros(length(window250_spikes),2);
timewindow_lfp = zeros(length(window250_spikes),2);

% Create a transition time window which can shift the window250_spikes to
% the corresponding time window in the lfp data.
for i= 1:length(window250_spikes)
    % Using abs() and min(), it finds the difference of the number that's 
    % closest to window250_spikes(1,i). 
    % Iterate through every column of the -250ms timestamps
    % And create a new variable (window250_lfp) containing these
    % timeshifts (in columns)
    window250_lfp(i,1) = min(abs(lfp_ts-window250_spikes(i,1)));
    window250_lfp(i,2) = min(abs(lfp_ts-window250_spikes(i,2)));
    
    % Create an lfp timewindow of -250ms to +250ms around spikes
    timewindow_lfp(i,1) = window250_spikes(i,1)+window250_lfp(i,1);
    timewindow_lfp(i,2) = window250_spikes(i,2)+window250_lfp(i,2);
end

% First for loop: iterates through the different timewindows (of -250 and
% 250 ms)
for j=1:length(timewindow_lfp) 
    % Find the positions corresponding to the values of the -250 and
    % 250 ms time window.
    tmin_pos = find(LFP_table.lfp_ts==timewindow_lfp(j,1));
    tmax_pos = find(LFP_table.lfp_ts==timewindow_lfp(j,2));
    %Preallocate Volt. This way, it is emptied during every iteration of j.
    Volt = [];
    % Second for loop: iterates through the values of lfp_ts and checks if
    % they fall within the j-th timewindow. 
    % Iteration from tmin_pos to tmax_pos.
    for i= tmin_pos:tmax_pos 
        if LFP_table.lfp_ts(i)>=timewindow_lfp(j,1) && LFP_table.lfp_ts(i)<=timewindow_lfp(j,2)
            % Create variable containing all Volt values within the time
            % window. 
            Volt(i,1)= LFP_table.lfp_data(lfp_ts==LFP_table.lfp_ts(i));
        end
    end
    % Create variable containing the means of each row 
    Volt_means(j) = mean(Volt(:,1));
end






% Stimulus presentation window
% is timeframe onset-offset of stimulus.




% Baseline window
% is always 500 ms before stimulus onset.



























