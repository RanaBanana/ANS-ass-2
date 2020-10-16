%Assignment 2
%Nick, Eline, & Rana 
%Nervous Neurons
%% Exercise 1
% Load assignment2_data.mat

% Combine the events_ts and events_type to mark when the stimulus was
% presented.
events_table = table(events_ts, events_type);

% Set of variables. Set to microseconds to match data.
timeBin = 25000; % time bin in microseconds
preStim = 500000; % time before stimulus onset in microseconds
postStim = 1000000; % time after stimulus offset in microseconds

% Extract values from table corresponding to the on- and offset
onTimes = events_table.events_ts(events_type==1);
offTimes = events_table.events_ts(events_type==31);

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
remove_zero = PSTH_time ~= 0;
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


%% Exercise 3
% Compute and plot the LFP power spectrum (using the pwelch function) during both
% baseline (500 ms before stimulus onset) and stimulus presentation (onset-offset)
% Take care in computing the power spectrum using parameters
% appropriate to the sampling frequency (1 kHz) and duration of the signal. Also compute and plot a
% relative power spectrum, defined as relative power change (per frequency bin) between
% baseline and presentation periods. Which (if any) LFP frequency is significantly
% modulated by stimulus presentation? What is the origin of the sharp peaks in the spectra?

clear all;
load('assignment2_data.mat');

%Variable definitions

%sampling frequency
Fs = 1000;
%Creating vectors with the timestamps during which the stimulus was
%presented (both on and off)
events = table(events_ts, events_type);
onTimes = events.events_ts(events_type==1);
offTimes = events.events_ts(events_type==31);

%LFP signal during stimulus presentation
%Create a matrix(preallocating)
lfp_stim = zeros(length(onTimes),500);

for i = 1:length(onTimes)
    [~,stim_begin_idx] = min(abs(lfp_ts - onTimes(i)));
    [~,stim_end_idx] = min(abs(lfp_ts - offTimes(i)));
    lfp_stim(i,(1:(stim_end_idx-stim_begin_idx+1))) = lfp_data(stim_begin_idx:stim_end_idx); 
end

%If this all doesn't work:
%Perhaps take the pwelch in the for-loop, put every P_welch value in an 
%array, and average it after the for-loop.

mean_lfp_stim = mean(lfp_stim,1);
%f = (-Fs/2, Fs/2, nfft)
%nfft_stim = length(mean_lfp_stim);
%f_stim = linspace(-Fs/2, Fs/2, nfft_stim);
[stim_P_welch, stim_F_welch] = pwelch(mean_lfp_stim,200,100,1000, Fs);
 
%LFP signal during baseline
%Create a matrix(preallocating)
lfp_bas = zeros(length(onTimes),500);

for j = 1:length(onTimes)
    [~,bas_begin_idx] = min(abs(lfp_ts - (onTimes(j)-500000)));
    [~,bas_end_idx] = min(abs(lfp_ts - onTimes(j)));
    lfp_bas(j,(1:(bas_end_idx-bas_begin_idx+1))) = lfp_data(bas_begin_idx:bas_end_idx); 
end

mean_lfp_bas = mean(lfp_bas,1);

[bas_P_welch, bas_F_welch] = pwelch(mean_lfp_bas,200,100,1000,Fs);

%Plotting stimulus presentation and baseline welch estimations
fig_lfp_stim = subplot(311);
plot(stim_F_welch,(log10(stim_P_welch)*10))
ylabel('Power');
xlim(fig_lfp_stim, [0 260])
xlabel('Normalised Frequency');
title("Welch Estimation during stimulus presentation")

fig_lfp_bas = subplot(312);
plot(bas_F_welch, (log10(bas_P_welch)*10))
ylabel('Power');
xlim(fig_lfp_bas, [0 260])
xlabel('Normalised Frequency');
title("Welch Estimation during baseline")

fig_lfp_bas_stim = subplot(313);
plot(stim_F_welch,(log10((stim_P_welch./bas_P_welch))*10))
ylabel('Power');
xlim(fig_lfp_bas_stim, [0 260])
xlabel('Normalised Frequency');
title("Welch Estimation for the relative power change (stimulus presentation/baseline)")


%% Exercise 4
% Exercise 4
% clear all
% load('assignment2_data.mat')

% Combine the events_ts and events_type to mark when the stimulus was
% presented.
events_table = table(events_ts, events_type);

% Extract values from table corresponding to the on- and offset 
onTimes = events_table.events_ts(events_type==1);
offTimes = events_table.events_ts(events_type==31);

% Set of variables. Set to microseconds to match data.
timeBin = 25000; % time bin in microseconds
ms_500 = 500000; % 500 ms in microseconds. Used as preStim and as stimulus presentation length. 
postStim = 1000000; % time after stimulus offset in microseconds
window_spike = 250000; % spike window -250 to 250ms around spikes

% Categorize the spikes in baseline and stimulus presentation window.
K=0;
L=0;
for i = 1:length(onTimes) 
    % 1. Baseline (-500 to 0 ms):
    % Find spikes_idx and corresponding ts in every baseline and add them
    % in one variable (base_spikes_ts)
    base_spikes_idx = find(spikes_ts>=(onTimes(i)-ms_500)& spikes_ts < onTimes(i));
    base_spikes_ts(K+1:K+numel(base_spikes_idx)) = spikes_ts(base_spikes_idx);
    K=K+numel(base_spikes_idx);
    
    % 2. Stimulus presentation window (0 to 500ms):
    % Find spikes_idx and corresponding ts in every stim presentation and 
    % add them in one variable (stim_spikes_ts) 
    stim_spikes_idx = find(spikes_ts>=onTimes(i)& spikes_ts < offTimes(i));
    stim_spikes_ts(L+1:L+numel(stim_spikes_idx)) = spikes_ts(stim_spikes_idx);
    L=L+numel(stim_spikes_idx);
end

% Find the corresponding lfp_ts values with the average spike window in baseline.
    % Preallocate window variables
    base_window_idx = NaN(length(base_spikes_ts), 2);
    base_window_length = NaN(length(base_spikes_ts),1);
    stim_window_idx = NaN(length(stim_spikes_ts), 2);
    stim_window_length =  NaN(length(stim_spikes_ts),1);    

% 1. Baseline:     
for j = 1: length(base_spikes_ts)
    [~, idx_beg] = min(abs(lfp_ts-(base_spikes_ts(j)-window_spike)));
    [~, idx_end] = min(abs(lfp_ts-(base_spikes_ts(j)+window_spike)));
    % Calculate the window length. This will later be used to calculate the
    % average time window. Idem for stimulus window.
    base_window_idx (j,:) = [idx_beg, idx_end]; 
    base_window_length(j) = idx_end-idx_beg;
end

% Calculate the average time window and iterate over the lfp_data signal to
% find the corresponding lfp values. The values start from the idx begin
% + the average length of the window-1. Idem for stimulus window.
base_window_length_mean_r = round(mean(base_window_length));
    
% Preallocate base_lfp
    base_lfp = NaN(length(base_window_idx), base_window_length_mean_r);
for k = 1:length(base_window_idx)
    base_lfp(k,:) = lfp_data(base_window_idx(k,1):(base_window_idx(k,1)+base_window_length_mean_r)-1);
end

% 2. Stimulus presentation (same steps as baseline):
for l = 1: length(stim_spikes_ts)
    [~, idx_beg] = min(abs(lfp_ts-(stim_spikes_ts(l)-window_spike)));
    [~, idx_end] = min(abs(lfp_ts-(stim_spikes_ts(l)+window_spike))); 
    stim_window_idx (l,:) = [idx_beg, idx_end]; 
    stim_window_length(l) = idx_end-idx_beg;
end

stim_window_length_mean_r = round(mean(stim_window_length));  
    
% Preallocate stim_lfp
    stim_lfp = NaN(length(stim_window_idx), stim_window_length_mean_r);
for m = 1:length(stim_window_idx)
    stim_lfp(m,:) = lfp_data(stim_window_idx(m,1):(stim_window_idx(m,1)+stim_window_length_mean_r)-1);
end

% Plot spike-triggered LFP signal for baseline and stimulus presentation:

% Plot lfp_base_mean:
figure('Name','Exercise 4','NumberTitle','off')
subplot(2,1,1)
lfp_base_mean = mean(base_lfp, 1);
t1 = linspace(-250,250, length(lfp_base_mean));
plot(t1,lfp_base_mean);
title('Spike-triggered LFP signal during baseline')
xlabel('Time window around spikes (ms)');
ylabel('Average LFP signal (µV)');

% Plot lfp_stim_mean:
subplot(2,1,2)
lfp_stim_mean = mean(stim_lfp, 1);
t2 = linspace(-250,250, length(lfp_stim_mean));
plot(t2,lfp_stim_mean);
title('Spike-triggered LFP signal during stimulus presentation')
xlabel('Time window around spikes (ms)');
ylabel('Average LFP signal (µV)');

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













