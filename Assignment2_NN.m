%Assignment 2
%Nick, Eline, & Rana 
%Nervous Neurons

%% Exercise 1
clear all

% load data
load('assignment2_data.mat')

% Combine the events_ts and events_type to mark when the stimulus was
% presented.
events_table = table(events_ts, events_type);

% Extract values from table corresponding to the on- and offset 
onTimes = events_table.events_ts(events_type==1);
offTimes = events_table.events_ts(events_type==31);

% Set of variables. Set to microseconds to match data.
ms_500 = 500000; % 500 ms in microseconds. Used as preStim and stimulus presentation length. 
postStim = 1000000; % time after stimulus offset in microseconds

% Select PSTH binsize and edges
timeBin = 25000; % time bin in microseconds
edges = (-500000: timeBin: 1500000);

% PSTH
spike_count=0;
for i = 1:length(onTimes) % Loop over all trials
    % Select spikes during the i-th trial
    hit = find(spikes_ts >= (onTimes(i) - ms_500) & spikes_ts < (offTimes(i) + postStim));
    
    % Calculate spike times relative to stimulus onset
    relative_spikes = spikes_ts(hit) - onTimes(i);
    
    % Calculate number of spikes per bin (25 ms). histcounts() divides the 
    % relative_spikes into bins specified by edges.
    spike_count = spike_count+histcounts(relative_spikes,edges);  
end

% Divide by total length of PSTH to compute firing rate (= estimated spikes 
% per second. In Hz, so timeBin in seconds)
firing_rate = spike_count/(length(onTimes)*(timeBin/10^6));

% Plot PSTH
figure('Name','Exercise 1','NumberTitle','off')
t_freq = linspace(-500,1500, length(firing_rate));
bar(t_freq,firing_rate);
xlabel('Time from stimulus onset (s)', 'FontSize', 12);
ylabel('Firing Rate (Hz)', 'FontSize', 12);
title('PSTH of audiovisual stimulus', 'FontSize', 14);
 
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
% Which (if any) LFP frequency is significantly
% modulated by stimulus presentation? What is the origin of the sharp peaks in the spectra?

%clearing workspace and loading raw data.
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
%Creating a matrix(preallocating) for the LFP data values during both
%stimulus presentation and baseline periods.
lfp_stim = zeros(length(onTimes),500);
lfp_bas = zeros(length(onTimes),500);

% This for-loop makes sure that the onset and offset timestamps correspond
% with the same LFP timestamps for the stimulus presentation periods. Also,
% it makes sure that the onset-500ms and onset timestamps correspond
% with the same LFP timestamps for the baseline periods. The iteration runs
% for length(onTimes) times, because this is the same number of stimulus
% presentation (and baseline) periods.
% Next, this for-loop takes the index for the wanted timestamps, which
% corresponds to the index in lfp_data. The wanted values in lfp_data are
% then cast into the matrix that was made before (lfp_stim or lfp_bas) for
% both the stimulus presentation and baseline periods.
for i = 1:length(onTimes)
    [~,stim_begin_idx] = min(abs(lfp_ts - onTimes(i)));
    [~,stim_end_idx] = min(abs(lfp_ts - offTimes(i)));
    lfp_stim(i,(1:(stim_end_idx-stim_begin_idx+1))) = lfp_data(stim_begin_idx:stim_end_idx);
    
    [~,bas_begin_idx] = min(abs(lfp_ts - (onTimes(i)-500000)));
    [~,bas_end_idx] = min(abs(lfp_ts - onTimes(i)));
    lfp_bas(i,(1:(bas_end_idx-bas_begin_idx+1))) = lfp_data(bas_begin_idx:bas_end_idx); 
end

% Calculating the mean LFP signal for the stimulus presentation period.
% Then, the pwelch function is used to make a Welch estimation.
mean_lfp_stim = mean(lfp_stim,1);
[stim_P_welch, stim_F_welch] = pwelch(mean_lfp_stim,200,100,1000, Fs);

% Calculating the mean LFP signal for the baseline period.
% Then, the pwelch function is used to make a Welch estimation.
mean_lfp_bas = mean(lfp_bas,1);
[bas_P_welch, bas_F_welch] = pwelch(mean_lfp_bas,200,100,1000,Fs);

% Plotting the Welch estimation for the stimulus presentation period.
fig_lfp_stim = subplot(311);
plot(stim_F_welch,(log10(stim_P_welch)*10))
ylabel('Power');
xlim(fig_lfp_stim, [0 260])
xlabel('Normalised Frequency');
title("Welch Estimation during stimulus presentation")

% Plotting the Welch estimation for the baseline period.
fig_lfp_bas = subplot(312);
plot(bas_F_welch, (log10(bas_P_welch)*10))
ylabel('Power');
xlim(fig_lfp_bas, [0 260])
xlabel('Normalised Frequency');
title("Welch Estimation during baseline")

% Plotting the Welch estimation for the relative power change (calculated
% by element-wise division of the power for stimulus presentation by the
% power for baseline.
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
% Rough expectations of what is desired.
% Consider the following 3 frequency ranges (2,6)Hz, (10, 20)Hz and (30,
% 40)Hz

% We need to evaluate: Spike Coherence in the ranges

% Compute if AP's happen more in certain phases in the LFP oscillation of
% that range of Hz.

%% 5a - Use an IIR filter to divide the LFP into these three frequency ranges.
% clear all; - remove comment notation if data is missing or mismatched :)
% close all; - Idem
% load assignment2_data.mat; - Idem

% Declare certain basic variables: frequency sampling, length of data etc
% etc
Fs=1000;
x = lfp_data; % So we have a simple var for our data
t = 0:length(x)-1; % we count from zero okay >=(

% We want to set up a few measures first that are relevant for a fourier 
% transform. Listed below are:
nfft = length(lfp_data); % Simple answer: Number of relevant datapoints
                         % Complex answer: Nfft is required to zero-pad our
                         % time-domain vector before we calculate any
                         % transformation (this makes it more efficiënt in
                         % the long run).
f=(-nfft/2:nfft/2-1)*Fs/nfft; % Specify our total range, multiplied by the
                              % fraction of sampling freq over nfft.
f_s= fftshift(f); % Neat way of shifting our spectrum's zero-frequency com-
                  % ponent to the center of our array.
X = fft(x, nfft); % Computes the -discrete Fourier transform- of our x values
                  % with the fast Fourier transform algorithm. Mathworks
                  % holds some interesting examples of how it works, if in-
                  % terested!

% We create relevant filters for the requested ranges from our assignment.
filt_2_6 = zeros(size(X)); % First we allocate space from our size of X.
a=find(abs(f_s)>= 2 & abs(f_s) <= 6); % We use the find() function again
                                      % to look for all points inbetween
                                      % the specified values.
filt_2_6(a) = 1; % Odd step at cursory glance, we assign 1 to "hits" so  
                 % that during later element operation we can pull a nifty
                 % trick (see line 57).

filt_10_20 = zeros(size(X)); % We repeat these steps twice below. 
b=find(abs(f_s)>= 10 & abs(f_s) <= 20); % -Same logic-
filt_10_20(b) = 1; % -Same logic- 

filt_30_40 = zeros(size(X)); % -Same logic-
c=find(abs(f_s)>= 30 & abs(f_s) <= 40); % -Same logic-
filt_30_40(c) = 1; % -Same logic-

X_f1 = X.*filt_2_6; % This was the nifty trick (Nick's imo, don't judge my
                    % group if it's lame). We want to multiply all relevant
                    % datapoints our filter with our Fourier transformed
                    % data: X. We use the .* application to multiply
                    % element wise through all data, because we assigned a
                    % 1 earlier, we only retain any points of relevance.
                    % (I suspect there is a cleaner method to this, as this
                    % does require an extra iteration of N calculations,
                    % but I couldn't come up with one in a reasonable
                    % timeframe.)
X_ff1 =ifft(X_f1, nfft); % We do an inverse fast Fourier transform to convert
                         % all data from our frequency domain to a time do-
                         % main. (if our signal is non-periodic, the resul-
                         % ting frequency spectrum will start to affected
                         % by leakage).
X_f2 = X.*filt_10_20; % Here we repeat the two steps above twice again.
X_ff2 =ifft(X_f2, nfft); % -Same logic-
X_f3 = X.*filt_30_40; % -Same logic-
X_ff3 =ifft(X_f3, nfft); % -Same logic-
Total_X_ff = {X_ff1 X_ff2 X_ff3}; % -Same logic

Freq = ["LFP 2 to 6Hz filtered", "LFP 10 to 20Hz filtered", "LFP 30 to 40Hz filtered"];
% ^ Assigning our desired title prints above for frequency is João's solu-
% tion, so credit goes to him! This variable will be used in the for-loop 
% below.

figure() ; % Call on a figure()
for plotId = 1 : 3 % We'll run it 3 times, as we desire 3 plots.
   subplot(3, 1, plotId) ; % Shift it with 1, per iteration (so we don't   
                           % create unnecessary overlap).
   plot(t(1:nfft), Total_X_ff{:, plotId}); % Plot for a range of 2000 points
                           % This value was chosen so the pattern is
                           % clearly shown to be repeating. 
   ylabel('mV');
   xlabel('datapoints');
   title(Freq(plotId)); % Here we call upon our Freq title variable with
                        % shifting position of (plotId).
end

%% 5b - Calculate the Hilbert transform for each filtered LFP trace

% We use the hilbert function to calculate the transform of our real input
% to a complext result of the same length. (useful for phases later)
y = hilbert(x); % We do this for the entire dataset, this is not asked 
                % by the exercise, but own interest demands it 8-).
y1 = hilbert(X_ff1); % We also do the transform for 2-6hz filter.
y2 = hilbert(X_ff2); % Same for the 10-20Hz variant.
y3 = hilbert(X_ff3); % Same for the 30-40Hz variant.

% We allocate a 1x4 cell for allowing a for-loop to iterate over the
% values of the hilbert transformation. 
y_tot = {y1 y2 y3 y};

% We declare another 1x4 string for easy title labeling later.
title_hil = ["Hilbert transform for 2 to 6Hz", "Hilbert transform for 10 to 20Hz", "Hilbert transform for 30 to 40Hz", "Hilbert transform for all Hz"];
    
figure() ; % Same method as used in exercise 5a, comments will be left out 
           % if the explanation is identical/similar enough as in 5a.
for plotId = 1 : 4
   subplot(4, 1, plotId) ;
   plot(t(1:2000), real(y_tot{:, plotId}(1:2000))); % We now plot over our 
                                                    % indexed cell struc-
                                                    % ture.
   hold on % Usage of hold on/off for overlapping the real-imag plots.
   plot(t(1:2000), imag(y_tot{:, plotId}(1:2000))); 
   hold off
   xlabel('Timeframe (ms)');
   ylabel('amplitude');
   title(title_hil(plotId)); % Same usage of iterating over the 1x4 string.
   legend('real element', 'imaginary element');
end


%% 5c - Calculate the instantaneous phase for the Hilbert transform.

% Declaring a 1x4 string again for later title usage.
title_phase = ["Hilbert phaseshift for 2 to 6Hz", "Hilbert phaseshift for 10 to 20Hz", "Hilbert phaseshift for 30 to 40Hz", "Hilbert phaseshift for original signal"];

% Same logic for this for-loop as the previous two. Only difference is 
% usage of angle() in the real-imag arguments. This way we get the rele-
% vant angle of all data relevant between specified filters.
figure() ;
for plotId = 1 : 4
   subplot(4, 1, plotId) ;
   plot(t(1:2000), real(angle(y_tot{:, plotId}(1:2000))));
   hold on
   plot(t(1:2000), imag(angle(y_tot{:, plotId}(1:2000))));
   hold off
   xlabel('Timeframe (ms)');
   ylabel('phase of angle');
   title(title_phase(plotId));
   legend('real element', 'imaginary element');
end

%% 5d & e - Per frequency range, plot a circular histogram of phases of spike time.

% Call upon three arguments below, as in exercise 1. Personally I just add
% these in case TA's review the exercises split seperately. 
events_table = table(events_ts, events_type);
onTimes = events_table.events_ts(events_type==1);
offTimes = events_table.events_ts(events_type==31);

% for-loop used to iterate over the length of stimuli. 
for k = 1:length(onTimes)
    % We assign any values relevant during stimulus presentation to "hit".
    hit = find(spikes_ts >= onTimes(k) & spikes_ts < offTimes(k));
    % Allocate to loc_sp any values of that position.
    loc_sp = spikes_ts(hit);
    % We create a m x n, of column (k) and rows of length(loc_sp).
    sp_dur_stim(k, 1:length(loc_sp)) = loc_sp;
end

remove_zero_sp = sp_dur_stim ~= 0; % remove any zero values left from our 
                                   % for-loop.
tryout = sp_dur_stim(remove_zero_sp)'; % Apply the removed zero's comparison
                                       % to our original m x n. Transpose
                                       % at the end for easier use later. 

% -Same logic as in 5a/b: Applies from line 174:200- Function calling would
% be useful in my opinion, but wasn't needed according to the teachers.
nfft2 = length(tryout); 
f1=(-nfft2/2:nfft2/2-1)*Fs/nfft2; 
f_s1= fftshift(f1); 

SDS = fft(tryout, nfft2); 
filt2_6 = zeros(size(SDS));  
a=find(abs(f_s1)>= 2 & abs(f_s1) <= 6);  
filt2_6(a) = 1; 

filt10_20 = zeros(size(SDS));
b=find(abs(f_s1)>= 10 & abs(f_s1) <= 20); 
filt10_20(b) = 1; 

filt30_40 = zeros(size(SDS)); 
c=find(abs(f_s1)>= 30 & abs(f_s1) <= 40); 
filt30_40(c) = 1; 

SDS2_6 = SDS.*filt2_6;
SDS1 = ifft(SDS2_6, nfft2);
SDS10_20 = SDS.*filt10_20;
SDS2 = ifft(SDS10_20, nfft2);
SDS30_40 = SDS.*filt30_40;
SDS3 = ifft(SDS30_40, nfft2);

ya1 = hilbert(SDS1);
ya2 = hilbert(SDS2);
ya3 = hilbert(SDS3);

% Simple plotting without the use for a for-loop, as the length of the con-
% struct doesn't get obscene here. Note the same usage of angle() as in 5c.
figure;
subplot(131);
polarhistogram(angle(ya1), 24);
title('Angle phases of spikes of 2 to 6Hz range during stimulus');
subplot(132);
polarhistogram(angle(ya2), 24);
title('Angle phases of spikes of 10 to 20Hz range during stimulus');
subplot(133);
polarhistogram(angle(ya3), 24);
title('Angle phases of spikes of 30 to 40Hz range during stimulus');
