%% Excercise 07-10-'20

f1 = 5; % Frequency 1 in 5 Hz
f2 = 50; % Frequency 2 in 50 Hz

%Sampling frequence has to conform to Nyquist theorem, so min of 2*20 =40.
Fs = 50; 
T_max = 10; 
nfft = 1024;
% Consider the following signals:
t = 0:1/Fs:T_max;

s1 = sin(2*pi*t*f1); % 5 Hz sine wave
s2 = sin(2*pi*t*f2); % 50 Hz sine wave
c = cos(2*pi*t*f2); % 5 Hz cosine wave

% % And their combinations:
% x1 = s1*0.8+s2*0.2;
% x2 = s1*0.8+c1*0.2;
% x3 = [s1,s2];

%plotting spectra of the waves
subplot(3,1,1)
plot(t,s1,t,s2,t,c)
xlabel('Time (s)')
title('Spectra of the waves')
legend('s1', 's2','c')

%plotting FFT magnitude spectra
f=linspace(-Fs/2,Fs/2,nfft);
S1=fft(s1,nfft);
S2=fft(s2,nfft);
C=fft(c,nfft);
subplot(3,1,2)
plot(f,abs(fftshift(S1)),f,abs(fftshift(S2)),f,abs(fftshift(C)))
xlabel('Frequency (Hz)')
title('FFT magnitude, shifted version')
legend('s1','s2','c')

%plotting FFT phase spectra
subplot(3,1,3)
plot(f,angle(fftshift(S1)),f,angle(fftshift(S2)),f,angle(fftshift(C)))
xlabel('Frequency (Hz)')
title('FFT phase, shifted version')
legend('s1','s2','c')

%% Results

%Spectra of the waves: s1 has a sine shape, whilst s2 and c are straight
% lines at 0 and 1 respectively. 
%FFT magnitude: s1 has a magnitude at peak 0, s2 is a horizontal line at 0, c
% has a magnitude peak at 0.
%Phase: s1 and c have similar phase between -5 and 5 Hz, s2 is slightly
% different than both.
% Fs = 50; 
% T_max = 0.1; 
% nfft = 1024;

%Spectra of the waves: s1 has a sine shape, whilst s2 and c are straight
% lines at 0 and 1 respectively. 
%FFT magnitude: s1 has a magnitutde peak at -5 and 5 Hz, s2 is a straight line,
% c has a magnitude peak at 0. There is a reasonably lot of leakage.
%Phase: s1 changes in amplitude between -5 and 5 Hz, s2 doesn't change, c
% goes 'lower' after 0 Hz.
% Fs = 50; 
% T_max = 2; 
% nfft = 1024;

%Spectra of the waves: s1 has a sine shape, with a higher frequency than
% previous settings. s2 and c are straight lines at 0 and 1 respectively. 
%FFT magnitude: s1 has two magnitude peaks at -5 and 5 Hz, s2 is a straight
% line (doesn't have a peak), c has a peak at 0 Hz. There is a lot less
% leakage than with the previous settings.
%Phase: s1 changes in amplitude between -5 and 5 Hz, s2 has a general 
% oscillating figure, c goes 'lower' after 0 Hz.
% Fs = 50; 
% T_max = 10; 
% nfft = 1024;

% Fs = 150; 
% T_max = 10; 
% nfft = 64;
%Spectra of the waves: All signals are shown as waves. However, detail of
% all waves is lost in 10 second timeframe.
%FFT magnitude: s1 has a magnitude around +5 and -5, s2 and s3 have magnitudes
% at -50 and +50 (both as expected). There is leakage (because of low nfft)
%Phase: s1 and s2 change in amplitude between -5 and +5. s2 and c change at 50Hz.

% Fs = 150; 
% T_max = 10;
% nfft = 1024;
%Spectra of the waves: All signals are shown as waves. However, detail of
% all waves is lost in 10 second timeframe.
%FFT magnitude: s1 has a magnitude around +5 and -5, s2 and s3 have magnitudes
% at -50 and +50 (both as expected). Little leakage because of long
% recording time.
%Phase: s1 and s2 change in amplitude between -5 and +5. s2 and c change at 50Hz.

% Fs = 1000; 
% T_max = 0.1; 
% nfft = 1024;
%Spectra of the waves: All signals are shown as waves (as a result of the
% high Fs)
%FFT magnitude: s1 has a magnitude around +5 and -5, s2 and s3 have magnitudes
% at -50 and +50 (both as expected). Leakage due to short recording
% time.
%Phase: s1 changes in amplitude between -5 and +5. 

%Spectra of the waves: All waves show sinuso?d behaviour or a variation of
% it. Important to note is that the frequency of  c and s2 are much higher
% that of s1, due to the frequency being at 50Hz instead of 5. 
%FFT magnitude: The magnitude of FFT is as expected, and relatively clean.
% however when zoomed in to their respective peaks the peaks appear wider.
% s1 shows peaks at -50/50, s2 is flat line at 0, c shows peaks at -5/5.
%Phase: s2 and c show an inverted wave of eachother but similar in other
% regards. s1 shows a different pattern and inverses itself between -5/5
% Fs = 1000; 
% T_max = 2; 
% nfft = 1024;

%Spectra of the waves: Incredibly dense spectrum is visible. Still follow 
% the same patterns as stated above. all run at a magnitude of -1/1
%FFT magnitude: s1 has a higher FFT magnitude in comparison to s2/c and
% peaks between 5 and 10Hz. s2/c peak at -40/55. 
%Phase: phase of s1 rises slowly when freq goes more towards positive, 
% however nearing 0 it rises more quickly. For s2 it follows the same slow 
% ascending pattern, only there is a single osscilating motion when
% reaching -50Hz. The opposite is true for c, where it follows a downward
% "escalator" motion from -50Hz.
% Fs = 1000; 
% T_max = 10; 
% nfft = 64;

%Spectra of the waves: No difference as the one directly above. As the they
% run the same Hz and Fs
%FFT magnitude: Amplitude of s2 and c peak both at -50/50. s1 peaks at -5/5
% Hz.
%Phase: c Shows a more drawn-out "escalator" motion from -50Hz to 50hz now.
% The same goes for s2, but now and inverse "escalator" motion. s1 is 
% relatively stable without much change until it reaches -5 to 5Hz, where 
% rapid ascension is gone through, before gradual descent, where rapid
% ascension happens at around 5Hz.
% Fs = 1000; 
% T_max = 10; 
% nfft = 1024;

%% Spectograms and filtering
% Set Fs=150; nfft=1024; T_max=10;
% Compute and display the spectrograms for x1,x2,x3 
t=0:1/Fs:2; % 2 secs @150Hz sample rate

% x1
y=chirp(t,5,1,50,'q'); % Start @ 5Hz, cross 50Hz at t=1sec
[s,F,T,P]=spectrogram(x1,128,120,128,1e3);

subplot(3,1,1)
surf(T,F,10*log(abs(P)), 'EdgeColor', 'none')
axis xy; axis tight; colormap(jet); view(0,90);
xlabel('Time (s)');
ylabel('Frequency(Hz)');

% x2
y=chirp(t,5,1,50,'q'); % Start @ 5Hz, cross 50Hz at t=1sec
[s,F,T,P]=spectrogram(x2,128,120,128,1e3);

subplot(3,1,2)
surf(T,F,10*log(abs(P)), 'EdgeColor', 'none')
axis xy; axis tight; colormap(jet); view(0,90);
xlabel('Time (s)');
ylabel('Frequency(Hz)');

% x3 = [s1,s2] goes from s1 with freq 5 to s2 with freq 50
y=chirp(t,5,1,50,'q'); % Start @ 5Hz, cross 50Hz at t=1sec
[s,F,T,P]=spectrogram(x3,128,120,128,1e3);% Computes the spectrogram of x3

subplot(3,1,3)

surf(T,F,10*log(abs(P)), 'EdgeColor', 'none')
axis xy; axis tight; colormap(jet); view(0,90);
xlabel('Time (s)');
ylabel('Frequency(Hz)');

%% Build a low pass filter in the frequency domain with a cutoff frequency of
%5.1 hz and one with 20Hz. 




