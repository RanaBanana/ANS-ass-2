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
PSTH = zeros(length(timeframePSTH));

for i = length(timeframePSTH) % idk man ik doe maar wat
PSTH = PSTH; % heel nuttig ook deze line
spikes

end

% Calculate the amount of spikes in each time bin. (prob not necessary?)
% N = histcounts(spikes_ts,timeBin);




%% Excercise 4
% Compute and plot the spike-triggered LFP signal (i.e. the average LFP [value of a] time signal
% occurring around the occurrence of an action potential) in the [-250, 250] ms window
% around the spike times. Do this separately for spikes in the baseline window and for
% spikes in the stimulus presentation window. What do you notice? Is there any difference
% in the LFP between baseline and stimulus presentation?

% Variable definitions
window_spike = 250000;
% lfp_data is expressed in uV.
LFP = [lfp_data ; lfp_ts];
%Test on first 12000 columns of LFP:
small_data = LFP(:,1:12000);
plot(small_data(2,:),small_data(1,:));

% Plotting the whole LFP signal; WARNING: ONLY RUN ON FAST COMPUTER
% plot_LFP = plot(LFP(2,:), LFP(1,:))

% Stimulus presentation window
% is timeframe onset-offset of stimulus.




% Baseline window
% is always 500 ms before stimulus onset.



























