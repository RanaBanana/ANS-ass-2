%% Assignment 2: Whole lotta data
% Anne Wijnen & Lois de Groot
clear all;
load('assignment2_data.mat');
%% Exercise 1:
onset_ts = events_ts(1:2:end);
offset_ts = events_ts(2:2:end);
begin_ts = onset_ts - 500*10^3;
ending_ts = offset_ts + 1000*10^3;

% subplot(211);
% hold on %Allow multiple plots on the same graph
% for j=1:length(spikes_ts) %Loop through each trial
%     t=spikes_ts(j);
%     for i = 1:length(t) %Loop through each spike time
%         line([t(i) t(i)], [j-1 j]) %Create a tick mark at x = t1(i) with a height of 1
%     end
% end
% xlabel('Time (sec)') %Label x-axis
% ylabel('Trials') %Label y-axis

%PSTH:
bin = 25*10^3;
edges = [-1:bin:1]; %Define edges of histogram
psth = zeros(length(edges),1); %Initialize the PSTH with zeros
for j=1:length(spikes_ts) %Loop over all trials
    %Add current trial?s spike times
    psth = psth + histc(spikes_ts(j),edges);
end
psth = psth/(length(spike)*bin_width);
subplot(212);
bar(edges,psth); %Plot PSTH as a bar graph
xlim([-1.1 1]) %Set limits of X-axis
xlabel('Time (sec)') %Label x-axis
ylabel('Firing rate (Hz)') %Label y-axis

% Responses
%Define baseline
bsl_win=find(edges<=-0.7);
%Compute basl mean
bsl_mean=mean(psth(bsl_win));
%Compute threshold for activation
bsl_std=std(psth(bsl_win));
th=bsl_mean+3*bsl_std;


%% Exercise 2:
onset_ts = events_ts(events_type == 1);
stim_ending = onset_ts + 500*10^3; %Zelfde als offset?
for i=1:length(spikes_ts)
    stim_idx = find(spikes_ts >= onset_ts(i) & spikes_ts < stim_ending(i));
    stim_win = spikes_ts(stim_idx);
    begin_ts = onset_ts(i) - 500*10^3;
    
    bas_idx = find(spikes_ts >= begin_ts & spikes_ts < onset_ts(i));
    bas_win = spikes_ts(bas_idx);
end
[h,p,ci,stats] = ttest(stim_win, bas_win)

%% Exercise 3:
nfft = length(lfp_data)*2
dur_signal = lfp_ts(end)-lfp_ts(1);
t = 0:1/lfp_fs:dur_signal; %Time vector of duration signal
LFP_data = fft(lfp_data,nfft);
[P_welch, F_welch]=pwelch(lfp_data,lfp_fs,[],nfft,lfp_fs);

%% Exercise 4:
clear all;
load('assignment2_data.mat');
figure;
plot(lfp_ts,lfp_data)


events = [events_type, events_ts];
onset_ts = events_ts(events_type == 1);
offset_ts = events_ts(events_type == 31);

spikes_ts_bas = [];
spikes_ts_stim = [];
spikes_baseline = [];
spikes_stim = [];
for i = 1:1000 %length(onset_ts)
    begin_ts = onset_ts(i) - 500000;
    baseline_idx = find(spikes_ts >= begin_ts & spikes_ts < onset_ts(i)); %wrm geen spikes_ts(i)
    spikes_ts_bas = spikes_ts(baseline_idx);
    spikes_baseline = horzcat(spikes_baseline, spikes_ts_bas);
    
    stim_idx = find(spikes_ts >= onset_ts(i) & spikes_ts < offset_ts(i));
    spikes_ts_stim = spikes_ts(stim_idx);
    spikes_stim = horzcat(spikes_stim, spikes_ts_stim); %wrm is dit nu wel veel waardes ipv spikes_stim?
end

%Length(spikes_ts)x2-matrix of the spike windows:
spikes_window1 = spikes_ts - 250*10^3;
spikes_window2 = spikes_ts + 250*10^3;
spikes_windows = [spikes_window1', spikes_window2'];

ERP_spikes_stim = [];
for i = 1:length(spikes_stim)
    [~, indexbeg] = min(abs(lfp_ts - (spikes_stim(i) - 250*10^3)));
    [~, indexend] = min(abs(lfp_ts - (spikes_stim(i) + 250*10^3)));
    ERP_spikes_stim(i,:) = lfp_data(indexbeg:indexend);
end

ERP_spikes_baseline = [];
for i = 1:length(spikes_baseline)
    [~, indexbeg] = min(abs(lfp_ts - (spikes_baseline(i) - 250*10^3)));
    [~, indexend] = min(abs(lfp_ts - (spikes_baseline(i) + 250*10^3)));
    ERP_spikes_baseline(i,:) = lfp_data(indexbeg:indexend);
end

% plotting ERP_baseline:
ERP_baseline = mean(ERP_spikes_baseline, 1);
plot(ERP_baseline)

% plotting ERP_stim:
ERP_stim = mean(ERP_spikes_stim, 1);
plot(ERP_stim)

%% Exercise 5:




