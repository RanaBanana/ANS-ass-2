% Assignment 2
% Nick, Eline, & Rana 
% Nervous Neurons

%% Exercise 1
clear all

% Load data
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
xlabel('Time from stimulus onset (ms)', 'FontSize', 12);
ylabel('Firing Rate (Hz)', 'FontSize', 12);
title('PSTH of audiovisual stimulus', 'FontSize', 14);
 
%% Exercise 2 --> continuation of exercise 1

% Create variable for 500000 ms.
ms_500 = 500000; % 500 ms in microseconds. Used as both time before and after stimulus onset.

% Categorize the firing rate bins in baseline and stimulus presentation 
% window.
baseline_fr = firing_rate(edges>=-ms_500 & edges<0); 
stimulus_fr = firing_rate(edges>=0 & edges<ms_500);

% Test assumptions for parametric testing
    % Test for normally distributed data using Shapiro-Wilk test:
    % H0: X is normal with unspecified mean and variance
[h_base,p_base,SW_base] = swtest(baseline_fr); % h = 1; p = 0.0240 <0.05; sw =0.8853-> reject H0
[h_stim,p_stim,SW_stim] = swtest(stimulus_fr); % h = 1; p = 0.0097 <0.05; sw =0.8653-> reject H0
    % Both are not normally distributed, so assumptions for parametric testing
    % are violated. --> non-parametric test: Wilcoxon singed rank test
    
% Wilcoxon Signed Rank Test:
    % H0 = x – y comes from a distribution with zero median.
[p_WSR, h_WSR, stats_WSR] = signrank(baseline_fr, stimulus_fr); % h = 1; p = 1.2975e-04 <0.05 --> reject H0

% This means that spiking activity during stimulus presentation period 
% significantly differs from spiking activity during baseline.

% Descriptive statistics:
    % Median
baseline_fr_med = median(baseline_fr); % median = 0.3640
stimulus_fr_med = median(stimulus_fr); % median = 14.7219

    % IQR
baseline_fr_IQR = iqr(baseline_fr); % IQR = 0.1011
stimulus_fr_IQR = iqr(stimulus_fr); % IQR = 14.2164

%% Exercise 3

% Clear workspace and loading raw data.
clear all;
load('assignment2_data.mat');

% Variable definitions
% Create vectors with the timestamps during which the stimulus was
% presented (for both on- and offset of stimulus presentation).
events = table(events_ts, events_type);
onTimes = events.events_ts(events_type==1);
offTimes = events.events_ts(events_type==31);

% Create a zeros-matrix (preallocating) for the wanted values of lfp_data 
stim_P = zeros(length(onTimes),501);
bas_P = zeros(length(onTimes),501);

% This for-loop makes sure that the onset and offset timestamps correspond
% with the same LFP timestamps for the stimulus presentation periods. Also,
% it makes sure that the onset-500ms and onset timestamps correspond
% with the same LFP timestamps for the baseline periods. The iteration runs
% for length(onTimes) times, because this is the same number of stimulus
% presentation (and baseline) periods.
% Next, this for-loop takes the index for the wanted timestamps, which
% corresponds to the index in lfp_data. The power and frequency values of 
% the wanted values in lfp_data are calculated using the pwelch() function.
% The power values are then cast into the matrix that was made
% before (stim_P or bas_P) for both the stimulus presentation and baseline 
% periods. stim_F and bas_F are the same for every iteration, so no matrix
% is needed for these variables.
for i = 1:length(onTimes)
    [~,stim_begin_idx] = min(abs(lfp_ts - onTimes(i)));
    [~,stim_end_idx] = min(abs(lfp_ts - offTimes(i)));
    [stim_P(i,:), stim_F] = pwelch(lfp_data(stim_begin_idx:stim_end_idx),200,100,1000,1000);
    
    [~,bas_begin_idx] = min(abs(lfp_ts - (onTimes(i)-500000)));
    [~,bas_end_idx] = min(abs(lfp_ts - onTimes(i))); 
    [bas_P(i,:), bas_F] = pwelch(lfp_data(bas_begin_idx:bas_end_idx),200,100,1000,1000);
end

% Calculating the mean power for the stimulus presentation period.
% Also, stim_F is assigned to a transposed version of stim_F, just for
% our peace of mind. This makes stim_F into a row vector instead of a
% column vector, which is nicer to work with.
mean_stim_P = mean(stim_P, 1);
stim_F = stim_F';

% Calculating the mean power for the baseline period.
% Also, bas_F is assigned to a transposed version of bas_F, just for
% our peace of mind. This makes bas_F into a row vector instead of a
% column vector, which is nicer to work with.
mean_bas_P = mean(bas_P, 1);
bas_F = bas_F';

% Plotting the power spectrum/ Welch estimation for both the stimulus 
% presentation period (red) and the baseline period (black).
fig_welch = subplot(211);
plot(stim_F,(log10(mean_stim_P)*10), 'r')
hold on
plot(bas_F, (log10(mean_bas_P)*10), 'k')
hold off
ylabel('Power (a.u)');
xlim(fig_welch, [0 260])
xlabel('Normalised Frequency (Hz)');
title("Power spectrum for stimulus presentation vs baseline")
legend('stimulus presentation', 'baseline')

% Plotting the power spectrum/ Welch estimation for the relative power 
% change (calculated by element-wise division of the power for stimulus 
% presentation by the power for baseline).
fig_relative_welch = subplot(212);
plot(stim_F,(log10((mean_stim_P./mean_bas_P))*10))
ylabel('Power (a.u)');
xlim(fig_relative_welch, [0 260])
xlabel('Normalised Frequency (Hz)');
title("Power spectrum for the relative power change (stimulus presentation/baseline)")

% Normalising the power vectors for baseline and stimulus presentation
normalise_baseline = mean_bas_P/ max(mean_bas_P);
normalise_stim = mean_stim_P/ max(mean_stim_P);

% Plotting the normalised power vs frequency
figure
plot(stim_F, 10*log10(normalise_baseline), 'k')
hold on
plot(stim_F, 10*log10(normalise_stim),'r')
ylabel('Normalised Power (a.u)');
xlim([0 260])
xlabel('Normalised Frequency (Hz)');
title("Power spectrum for stimulus presentation vs baseline")
legend('baseline', 'stimulus presentation')

% Statistics
% By visual inspection, we chose range1 to be from 10-80 Hz.
range1_stim = mean_stim_P(stim_F>10 & stim_F<=80);
range1_bas = mean_bas_P(stim_F>10 & stim_F<=80);

% By visual inspection, we chose range2 to be from 100-245 Hz.
range2_stim = mean_stim_P(stim_F>100 & stim_F<=245);
range2_bas = mean_bas_P(stim_F>100 & stim_F<=245);

% Assumptions
% Testing if the data is normally distributed
[H, pValue, W] = swtest(range1_stim);
% H = 1
% pValue = 2.1242e-09
% W = 0.6924

[H, pValue, W] = swtest(range1_bas);
% H =1
% pValue = 2.2443e-11
% W =0.5551

[H, pValue, W] = swtest(range2_bas);
% H =1
% pValue = 4.2839e-10
% W = 0.8332

[H, pValue, W] = swtest(range2_stim);
% H =1
% pValue =4.4787e-06
% W =0.9370

% The assumptions of normality of the data has been violated, so the 
% Wilcoxon signed rank test will be performed on our data.

% Wilcoxon signed rank test on range1
[p_WSR, h_WSR, stats_WSR] = signrank(range1_stim, range1_bas);
% p_WSR = 3.5595e-13 (p_value)
% h_WSR = 1
% stats_WSR = 
%   struct with fields:
%           zval: 7.2713
%     signedrank: 2485

% Wilcoxon signed rank test on range2
[p_WSR, h_WSR, stats_WSR] = signrank(range2_stim, range2_bas);
% p_WSR = 1.5247e-25 (p_value)
% h_WSR = 1
% stats_WSR = 
%   struct with fields:
%           zval: 10.4462
%     signedrank: 10585


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

%% 5 - Rough expectations of what is desired.
% Consider the following 3 frequency ranges (2,6)Hz, (10, 20)Hz and (30,
% 40)Hz

% We need to evaluate: Spike Coherence in the ranges

% Compute if AP's happen more in certain phases in the LFP oscillation of
% that range of Hz.

%% 5a - Use an IIR filter to divide the LFP into these three frequency ranges.
clear all;
close all;
load assignment2_data.mat;

% Declare certain basic variables: frequency sampling, length of data etc
% etc
Fs=1000;
x = lfp_data; % So we have a simple var for our data
t = 0:length(x)-1; % count from zero 

% We want to set up a few measures first that are relevant for a fourier 
% transform. Listed below are:
nfft = length(lfp_data); % Simple answer: Number of relevant datapoints
                         % Complex answer: Nfft is required to zero-pad our
                         % time-domain vector before we calculate any
                         % transformation (this makes it more efficient in
                         % the long run).
f=(-nfft/2:nfft/2-1)*Fs/nfft; % Specify our total range, multiplied by the
                              % fraction of sampling freq over nfft.
f_s= fftshift(f); % Neat way of shifting our spectrum's zero-frequency com-
                  % ponent to the center of our array.
X = fft(x, nfft); % Computes the -discrete Fourier transform- of our x values
                  % with the fast Fourier transform algorithm. Mathworks
                  % holds some interesting examples of how it works, if in-
                  % terested!

% Create relevant filters for the requested ranges from the assignment.
filt_2_6 = zeros(size(X)); % First we allocate space from our size of X.
a=find(abs(f_s)>= 2 & abs(f_s) <= 6); % We use the find() function again
                                      % to look for all points inbetween
                                      % the specified values.
filt_2_6(a) = 1; % Odd step at cursory glance, we assign 1 to "hits" (a) so  
                 % that during later element operation we can pull a nifty
                 % trick (see line 388).

filt_10_20 = zeros(size(X)); % We repeat these steps twice below. 
b=find(abs(f_s)>= 10 & abs(f_s) <= 20);
filt_10_20(b) = 1; 

filt_30_40 = zeros(size(X)); 
c=find(abs(f_s)>= 30 & abs(f_s) <= 40); 
filt_30_40(c) = 1; 

X_f1 = X.*filt_2_6; % This was the nifty trick. We want to multiply all relevant
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
                         % main. (If the signal is non-periodic, the resul-
                         % ting frequency spectrum will start to affected
                         % by leakage).
X_f2 = X.*filt_10_20; % Here we repeat the two steps above twice again.
X_ff2 =ifft(X_f2, nfft); 
X_f3 = X.*filt_30_40; 
X_ff3 =ifft(X_f3, nfft); 
Total_X_ff = {X_ff1 X_ff2 X_ff3};

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
                % by the exercise, but own interest demands it.
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

% -Same logic as in 5a/b: Applies to line 507:538- Function calling would
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

% Use remove zero logic, and iterate before transposing. 
% Logic applies same way as in 5a/b. 
SDS2_6 = SDS.*filt2_6;
SDS2_6_rem0 = SDS2_6 ~= 0;
SDS2_6 = SDS2_6(SDS2_6_rem0)';
SDS1 = ifft(SDS2_6, nfft2);
SDS10_20 = SDS.*filt10_20;
SDS10_20_rem0 = SDS10_20 ~= 0;
SDS10_20 = SDS10_20(SDS10_20_rem0)';
SDS2 = ifft(SDS10_20, nfft2);
SDS30_40 = SDS.*filt30_40;
SDS30_40_rem0 = SDS30_40 ~= 0;
SDS30_40 = SDS30_40(SDS30_40_rem0)';
SDS3 = ifft(SDS30_40, nfft2);

% Hilbert transform of SDS results, after that angle() is applied so we
% use a different variable name (in angles) for the circular statistics. 
ya1 = hilbert(SDS1);
ya2 = hilbert(SDS2);
ya3 = hilbert(SDS3);
ya_a1 = angle(ya1); 
ya_a2 = angle(ya2);
ya_a3 = angle(ya3);

circ_rtest(ya_a1); % p = 0
circ_rtest(ya_a2); % p = 1.636*10^-88 -> (0) :)
circ_rtest(ya_a3); % p = 0

cry1 = circ_kuipertest(ya_a1, ya_a2); % p = 0.001
cry2 = circ_kuipertest(ya_a1, ya_a3); % p = 0.001
cry3 = circ_kuipertest(ya_a2, ya_a3); % p = 0.001

% Simple plotting without the use for a for-loop, as the length of the con-
% struct doesn't get obscene here. Note the same usage of angle() as in 5c.
figure;
subplot(131);
polarhistogram(angle(ya1), 24);
title('Angle phases of 2 to 6Hz range during stimulus');
DisplayStyle = 'stairs'
subplot(132);
polarhistogram(angle(ya2), 24);
title('Angle phases of 10 to 20Hz range during stimulus');
subplot(133);
polarhistogram(angle(ya3), 24);
title('Angle phases of 30 to 40Hz range during stimulus');
DisplayStyle = 'stairs'
