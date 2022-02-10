clc;clear;close all;
% Read the two signals (messages) %
[m1,Fs1]=audioread('Short_QuranPalestine.wav');
[m2,Fs2]=audioread('Short_FM9090.wav');
% Transferring from two channels signal to one channel %
m1 = m1(:,1) + m1(:,2);
m2 = m2(:,1) + m2(:,2);
% Padding the short signal to be the same length as long signal %
m2 = [m2;zeros(length(m1) - length(m2),1)];
% plot fft of two messages on baseband (frequency domain) and check their BW %
M1=fft(m1);
k = -length(M1)/2 : length(M1)/2-1;
subplot(3,2,1)
plot (k*Fs1/length(M1),abs(fftshift(M1)))
title('Baseband spectrum of message 1')
xlabel('Frequency in Hz')
M2=fft(m2);
k = -length(M2)/2 : length(M2)/2-1;
subplot(3,2,2)
plot (k*Fs2/length(M2),abs(fftshift(M2)))
title('Baseband spectrum of message 2')
xlabel('Frequency in Hz')
% increasing sample frequency of each signal (Nyquist criteria) as Fc = 100KHz %
m1 = interp(m1,10);
Fs1 = Fs1 * 10;
m2 = interp(m2,10);
Fs2 = Fs2 * 10;
% AM modulator % 
Fc1=100*10^3; %frequency of first carrier
Fc2=150*10^3; %frequency of second carrier
T1 = 0:(1/Fs1): length(m1)/(Fs1) - 1/(Fs1);
T2 = 0:(1/Fs2): length(m2)/(Fs2) - 1/(Fs2);
carrier1 = cos((2*pi*Fc1)*T1)';
carrier2 = cos((2*pi*Fc2)*T2)';
m1 = m1.*carrier1;
m2 = m2.*carrier2;
% Plot the modulation of the signals in FD %
M1=fft(m1);
k = -length(M1)/2 : length(M1)/2-1;
subplot(3,2,3)
plot (k*Fs1/length(M1),abs(fftshift(M1)))
title('Modulation of message 1 with carrier frequrncy 100 KHz')
xlabel('Frequency in Hz')
M2=fft(m2);
k = -length(M2)/2 : length(M2)/2-1;
subplot(3,2,4)
plot (k*Fs2/length(M2),abs(fftshift(M2)))
title('Modulation of message 2 with carrier frequrncy 150 KHz')
xlabel('Frequency in Hz')
% Plot the spectrum which will be sent %
subplot(3, 2, [5 6]);
mx = m1 + m2;
Mx = fft(mx);
plot(k*Fs1/length(Mx), abs(fftshift(Mx)))
title('The spectrum of the output of the transmitter')
xlabel('Frequency in Hz')
%RF stage% 
% Preparing the BandPass Filter to each signal %
BPF1 = designfilt('bandpassfir', 'FilterOrder', 20, ...
             'CutoffFrequency1', 90000, 'CutoffFrequency2', 110000,...
             'SampleRate', 441000);  
BPF2 = designfilt('bandpassfir', 'FilterOrder', 20, ...
             'CutoffFrequency1', 140000, 'CutoffFrequency2', 160000,...
             'SampleRate', 441000); 
% Plot the signals after the BPF in FD %
figure
subplot(2, 1, 1)
mx_RF1 = fftfilt(BPF1, mx);
%mx_RF1 = mx;
Mx_RF_BPF1 = fft(mx_RF1);
plot(k*Fs1/length(Mx), abs(fftshift(Mx_RF_BPF1)))
title('Message 1 after the RF filter (before the mixer)')
xlabel('Frequency in Hz')
subplot(2, 1, 2)
mx_RF2 = fftfilt(BPF2, mx);
%mx_RF2 = mx;
Mx_RF_BPF2 = fft(mx_RF2);
plot(k*Fs1/length(Mx), abs(fftshift(Mx_RF_BPF2)))
title('Message 2 after the RF filter (before the mixer)')
xlabel('Frequency in Hz')
% Mixer and  Oscillator %
T1 = 0 : 1/Fs1 : length(m1)/Fs1 - 1/Fs1;
carrier1 = cos(2*pi*(Fc1+25000)*T1)';
carrier2 = cos(2*pi*(Fc2+25000)*T1)';
mx_RF1_IF = mx_RF1 .* carrier1;
mx_RF2_IF = mx_RF2 .* carrier2;
Mx_RF1_IF = fft(mx_RF1_IF);
Mx_RF2_IF = fft(mx_RF2_IF);
% Ploting the signals after the mixer %
figure
subplot(2, 1, 1)
plot(k*Fs1/length(Mx_RF1_IF), abs(fftshift(Mx_RF1_IF)))
title('Message 1 after the Mixer')
xlabel('Frequency in Hz')
subplot(2, 1, 2)
plot(k*Fs2/length(Mx_RF2_IF), abs(fftshift(Mx_RF2_IF)))
title('Message 2 after the Mixer')
xlabel('Frequency in Hz')
% IF Stage %
% Preparing the BPF to Pass the signal at low frequency IF %
BPF3_IF = designfilt('bandpassfir', 'FilterOrder', 30, ...
             'CutoffFrequency1', 15000, 'CutoffFrequency2', 35000,...
             'SampleRate', 441000); 
% Ploting the signal part at IF frequency %
figure
subplot(2, 1, 1)
mx_IF1 = fftfilt(BPF3_IF, mx_RF1_IF);
Mx_IF1 = fft(mx_IF1);
plot(k*Fs1/length(Mx_IF1), abs(fftshift(Mx_IF1)))
title('Message 1 after the IF filter ')
xlabel('Frequency in Hz')
subplot(2, 1, 2)
mx_IF2 = fftfilt(BPF3_IF, mx_RF2_IF);
Mx_IF2 = fft(mx_IF2);
plot(k*Fs1/length(Mx_IF2), abs(fftshift(Mx_IF2)))
title('Message 2 after the IF filter ')
xlabel('Frequency in Hz')
% Baseband detection stage %
% creating the base carrier one with freq = 30 khz %
T = 0 : 1/(Fs1) : length(m1)/(Fs1) - 1/(Fs1);  
c = cos(2*pi*(25000)*T)';
% Ploting the signal before the last LPF %
figure
mx_base1 = mx_IF1 .* c;
Mx_base1 = fft(mx_base1);
subplot(2, 1, 1);
plot(k*Fs1/length(Mx_base1), abs(fftshift(Mx_base1)))
title('Message 1 after the mixer (before the LPF)')
xlabel('Frequency in Hz')
mx_base2 = mx_IF2 .* c;
Mx_base2 = fft(mx_base2);
subplot(2, 1, 2);
plot(k*Fs2/length(Mx_base2), abs(fftshift(Mx_base2)))
title('Message 2 after the mixer (before the LPF)')
xlabel('Frequency in Hz')
% Preparing the LowPass Filter to recover the original signal %
LPF = designfilt('lowpassfir', 'PassbandFrequency', 0.05,...
             'StopbandFrequency', 0.07, 'PassbandRipple', 0.5, ...
             'StopbandAttenuation', 65, 'DesignMethod', 'kaiserwin'); 
% Ploting the signals at its baseband in the receiver %
figure     
subplot(2, 1, 1);         
m_base1 = fftfilt(LPF, mx_base1);   
M_base1 = fft(m_base1);
plot(k*Fs1/length(M_base1), abs(fftshift(M_base1))) 
title('Message 1 after the LPF')
xlabel('Frequency in Hz')
subplot(2, 1, 2);
m_base2 = fftfilt(LPF, mx_base2);   
M_base2 = fft(m_base2);
plot(k*Fs1/length(M_base2), abs(fftshift(M_base2)))
title('Message 2 after the LPF')
xlabel('Frequency in Hz')
% DownSampling and Playing the audios %
output_sig1 = downsample(m_base1, 10);
output_sig2 = downsample(m_base2, 10);
sound(10 * output_sig1, Fs1/10)
pause(17)
sound(10 * output_sig2, Fs2/10)