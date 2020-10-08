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

% Iterate over length of all onset stimulus points (offset could also be
% done but less intuitive)
for i = 1:length(onTimes) 
    % We use find() to acquire all relevant timepoints on the dataset that are inbetween
    % the onset timing - 500 and offset timing + 1000. We then assign it to
    % a variable called hit.
    hit = find(spikes_ts >= (onTimes(i) - preStim) & spikes_ts <= (onTimes(i) + postStim));
    % we assign a temporary value for the spikestiming - onsettiming 
    % (REMOVE THIS COMMENT AT THE END - we do this so we can get a general
    % timepoint each spike is measured *relative!* to the onset of the
    % stimulus, it gives a handy way to plot it later I hope)
    spikes_on_time = spikes_ts(hit) - onTimes(i);
    % has all timepoints per measures spike minus corresponding onset of
    % stimulus. The value that spikes_on_time get is also used to show the
    % length needed so we dont have to create alternative loops
    % (REMOVE THIS COMMENT AT THE END - HOLY SHIT THIS IS UGLY. It works
    % but we assign a big fucking chunky disaster of a method, as PSTH_time
    % changes every fucking time we have a different spikes/i so it's omega
    % inefficiënt. But I'm too braindead to fix this now, remind me when we
    % have time later or maybe one of you sees a good solution, for now
    % bear with this integrity violating construct pl0x 8-)).
    PSTH_time(i, 1:length(spikes_on_time)) = spikes_on_time;
end

% Easy way to change our microseconds to milliseconds.
PSTH_time = PSTH_time / 1000;

% We do this so we can get rid of the zeros, otherwise it fucks our whole
% graph because it will count the zeros as hits (this is a result of the
% moral corruptness mentioned above, it assigns 0's when it's filling in
% unneccesary space :'). so we fix it with a bandaid as shown below.
remove_zero = PSTH_time > 0;
% IMPORTANT NOTE: yo Nick your bandaid fucking sucks. Ik zag net op
% whatsapp dat de frequency count 300 verschilt met en zonder filter. Not
% sure how maar maybe ramt hij ontzettend kleine getallen bij 0 er ook uit.
% I'll take a look at it tomorrow

% Dumb name but we use to create a vector of 50 to 100 with so many fucking
% steps it gets high enough in our histogram to make the orange bar. I
% wanted it in and got stubborn to show where the stimulus duration was
% present. It's ugly now but nothing crying cant fix
testing = (50:0.08:100);

% All the outward stuff to make the histogram, first half stolen off
% Eline's code. Rest is self explanatory except for hold on/off maybe. We
% use it to be able to push two plots on top of eachother, they don't like
% eachother I'm pretty sure but using the hold on cheat code gave them no
% choice. If they start fighting blame Matlab please.
figure
histogram(PSTH_time(remove_zero), 80)
xlabel('Time (milliseconds)', 'FontSize', 12);
ylabel('Spike Count', 'FontSize', 12);
axis tight
title('Spike frequency per 25 ms time bin', 'FontSize', 14);
hold on
histogram(testing, 1)
hold off
% To be honest not sure why the legend function works as this. It's
% strangely fucking neat it assigns the appropriate colors to the calls.
% Kinda cool
legend('Spikes per timebin', 'timeframe of stimulus');



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



























