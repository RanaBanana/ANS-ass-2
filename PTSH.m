%Assignment 2
%% Excercise 1
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

%% Exercise 2
% Evaluate whether spiking activity during stimulus presentation period (the period
% from stimulus onset to 500 ms afterwards) is significantly different from baseline
% (defined as the activity in the 500 ms preceding stimulus onset). Is there a significant
% response?

%These are both the same value, so to make the code more efficient it is
%possible to make these the same variable.
preStim = 500000; % time before stimulus onset in microseconds
postStim = 500000; %time after stimulus onset in microseconds

for i = 1:length(onTimes) 
    % We use find() to acquire all relevant timepoints on the dataset that are inbetween
    % the onset timing - 500 and onset timing + 500. We then assign it to
    % a variable called hit.
    hit = find(spikes_ts >= (onTimes(i) - preStim) & spikes_ts <= (onTimes(i) + postStim));
    %Baseline
    hit_baseline = find(spikes_ts <= onTimes(i) & spikes_ts >= (onTimes(i) - preStim));
    %stimulus presentation
    hit_stm_pres = find(spikes_ts >= onTimes(i) & spikes_ts <= (onTimes(i) + postStim));
    
    % Create vectors containing the spiketimes in all baseline (-500 to
    % 0ms) and all stimulus (0 to 500) timewindow
    baseline(i, 1:length(spikes_on_time_bas)) = spikes_on_time_bas;
    stimulus(i, 1:length(spikes_on_time_stm)) = spikes_on_time_stm;
    
end

% Convert timestamps in both timewindows to 1s and 0s in order to count the
% different spikes. (This does not remove the 0s!)
baseline_remove_zeros = baseline ~=0;
stimulus_remove_zeros = stimulus ~=0;

% Calculate the #spikes for each timewindow in all baseline and stimulus
% timewindow. 
% !!!!!! This can probably be done in a different (quicker way), same
% for the step above.
baseline_small_sum = sum(baseline_remove_zeros,2);
stimulus_small_sum = sum(stimulus_remove_zeros,2);

% Calculate the average #spikes in the baseline and in the stimulus
% timewindow.
baseline_avg = mean(baseline_small_sum); % 0.1749
stimulus_avg = mean(stimulus_small_sum); % 6.0566


%% Excercise 3
% Compute and plot the LFP power spectrum (using the pwelch function) during both
% baseline (500 ms before stimulus onset) and stimulus presentation (onset-offset)
% Take care in computing the power spectrum using parameters
% appropriate to the sampling frequency and duration of the signal. Also compute a plot a
% relative power spectrum, defined as relative power change (per frequency bin) between
% baseline and presentation periods. Which (if any) LFP frequency is significantly
% modulated by stimulus presentation? What is the origin of the sharp peaks in the spectra?

LFP = [lfp_data ; lfp_ts];
min(abs(lfp_ts- events_ts))
%Stimulus presentation


%Baseline




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

%% exercise 5

%lmao fml
clear all;
close all;
load assignment2_data.mat;

% randomshit gestolen van slides/code/google/wanhoop
Fs=1000;
N = 20;
x = lfp_data; %deze nog wel zelfbedacht :')
t = 0:length(x)-1; %deze ook want ik ben te fragiel om 1:length(x) te doen.

% Compute fourier transform
nfft = length(lfp_data); % Ab
f=(-nfft/2:nfft/2-1)*Fs/nfft; %ra 
f_s= fftshift(f); %ca
X = fft(x, nfft); %da
%bruh

%Hier toen we 3x een filter voor de gewenste ranges
filt_2_6 = zeros(size(X)); % allocaten ruimte voor size van X
a=find(abs(f_s)>= 2 & abs(f_s) <= 6); % zoeken naar alle relevante waarden van f_s (freq)
filt_2_6(a) = 1; % if present we adjust the values to 1 for later (dan heb je dus een hit van die freq range)

filt_10_20 = zeros(size(X)); % idem
b=find(abs(f_s)>= 10 & abs(f_s) <= 20); %diedum
filt_10_20(b) = 1; %riedum

filt_30_40 = zeros(size(X)); %smiedum
c=find(abs(f_s)>= 30 & abs(f_s) <= 40); %geenrijmmeerdum
filt_30_40(c) = 1; % tsja, de comments zijn mager

X_f1 = X.*filt_2_6; %magical shit van onze filter voorwaarden(?) met onze x punten
X_ff1 =ifft(X_f1, nfft); %inverse fast fourier transform (ik pretendeer later wel te weten wat het is) over de 2 waarden
X_f2 = X.*filt_10_20;
X_ff2 =ifft(X_f2, nfft);
X_f3 = X.*filt_30_40;
X_ff3 =ifft(X_f3, nfft);

%plot tijd en kijken of poging nfft^length(lfp_data) wel werkt :') 
figure;
figfor2_6 = subplot(311);
plot(t(1:nfft), X_ff1);
ylim(figfor2_6, [-5000 5000]);
ylabel('mV');
xlabel('datapoints');
title('LFP 2 to 6 Hz filtered');

figfor10_20 = subplot(312);
plot(t(1:nfft), X_ff2);
ylim(figfor10_20, [-5000 5000]);
ylabel('mV');
xlabel('datapoints');
title('LFP 10 to 20 Hz filtered');

figfor30_40 = subplot(313);
plot(t(1:nfft), X_ff3);
ylim(figfor30_40, [-2000 2000]);
ylabel('mV');
xlabel('datapoints');
title('LFP 30 to 40 Hz filtered');
% courtesy of Rana https://gyazo.com/806acc8c7477862103ae523c9845507b



% Just kidding it was okay, it just feel wrong to have done it with just
% restructering copied code from slides but hey. 

























