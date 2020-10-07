%% Excercise 07-10-'20

f1 = 5; % Frequency 1 in 5 Hz
f2 = 50; % Frequency 2 in 50 Hz

%Sampling frequence has to conform to Nyquist theorem, so min of 2*20 =40.
Fs = 50; 
T_max = 0.1; 
nfft = 1024;
% Consider the following signals:
t = 0:1/Fs:T_max;

s1 = sin(2*pi*t*f1); % 5 Hz sine wave
s2 = sin(2*pi*t*f2); % 50 Hz sine wave
c = cos(2*pi*t*f2); % 5 Hz sine wave

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
% lines. 
%FFT magnitude: s1 has a magnitude at 0, s2 is a horizontal line at 0, c
% has a magnitude at 0.
%Phase: 
% Fs = 50; 
% T_max = 0.1; 
% nfft = 1024;






