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
