% 
% Práctica 2.2: "Signal Analysis: ASK simulation in Matlab with media effects."
% Rodo Vera
% •
% 17:39
% 100 puntos
% Fecha de entrega: 3 may, 16:00
% Práctica 2.2: "Signal Analysis: ASK simulation in Matlab with media effects."
% 
% PRACTICE 2.2.pptx
% PowerPoint
% 
% test1.m
% Objective C
% Comentarios de la clase
% Tu trabajo
% Asignado
% Comentarios privados
clear all;
clf;

load filterLP3.mat;

A = 20; %W transmitted power (43dBm) EXAMPLE: CELLULAR TOWER
fc = 10e3; %kHz; %Hz Carrier Frec

fs = 44100; %samples per second
Ts = 1/fs; %s

tmax = 2; %s 2 sec
t1 = 0:Ts:tmax; %0->tmaxs [s]
s1 = A*sin(2*pi*fc*t1);


figure(1);
plot(t1,s1);
xlabel('Time [s]');
title("Sine signal with f_c= " + fc + " [Hz]");
phase = pi;
s2 = A*sin(2*pi*fc*t1 + phase);

figure(2);
plot(t1,s2,'r');
xlabel('Time [s]');
title("Sine signal with a de-phase of 180°");
%%in the air:
%lambda = c / fc;


s3 = A/3*sin(2*pi*3*fc*t1);
s4 = A/5*sin(2*pi*5*fc*t1);

figure(3);
plot(t1,s3);
xlabel('Time [s]');
title("Sine signal with 3*f_c");

sTotal = s1 + s3 + s4;
figure(4);
plot(t1,sTotal);
xlabel('Time [s]');
title("Sine signal s1 + s3 + s4");

%%%Frequency analysis
N = fs*tmax; 
freq = (-fs/2):(fs/N):(fs/2);

%%%Frequency analysis Stotal
Stotal = abs(fftshift(fft(sTotal))); %DEP (PSD)
figure(5);
plot(freq,Stotal);
xlabel('Freq [Hz]');
title("Frequency Spectrum of sine signal s1+s2+s3");



figure(6);
S1 = abs(fftshift(fft(s1))); %DEP (PSD)
plot(freq,S1);
xlabel('Freq [Hz]');
title("Frequency Spectrum of sine signal with f_c= " + fc + " [Hz]");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% %%%MODULATING IN ASK SIMULATION:
pulsoBitUno(1:(length(t1)/2)+1)=A; %bit 1;
tpulso = t1;
pulsoBitCero(1:length(t1)/2)= 0; %Bit 0;


pulsos = horzcat(pulsoBitUno,pulsoBitCero);
figure(7);
plot(tpulso,pulsos);
xlabel('Time [s]');
title("Bit Pulses 1 and 0");



%%
modulatedSignalASK = s1 .* pulsos;


figure(8);
plot(tpulso,modulatedSignalASK);
xlabel('Time [s]');
title("ASK bit modulation with f_c= " + fc + "[Hz]");
%%%Frequency analysis of Modulated Signal
N = fs*tmax;
PULSOSfrec = abs(fftshift(fft(pulsos))); %DEP (PSD)


freq = (-fs/2):(fs/N):(fs/2);
figure(9);
plot(freq,PULSOSfrec);
xlabel('Freq [Hz]');
title("Frequency Spectrum of 1 and 0 bits");
modulatedSignalASKfrec = abs(fftshift(fft(modulatedSignalASK))); %DEP (PSD)

freq = (-fs/2):(fs/N):(fs/2);
figure(10);
plot(freq,modulatedSignalASKfrec,'r');
xlabel('Freq [Hz]');
title("Frequency Spectrum of ASK bit modulation with f_c= " + fc + "[Hz]");


%%%ADDING NOISE AWGN
SNR=25; %%dBs
modulatedSignalASKnoise = awgn(modulatedSignalASK,SNR,'measured');
figure(11);
plot(tpulso,modulatedSignalASKnoise);
xlabel('Time [s]');
title("ASK signal with awgn noise @" + SNR +"dB");


%%%ADDING INTERFERENCE (OTHER SIGNALS ADDING TO NOISE)
A1 = 20; A2 = 20; A3 = 20;
fc1 = 10000; fc2 = 900; fc3 = 500;

int1 = A1*cos(2*pi*fc1*t1);
int2 = A2*cos(2*pi*fc2*t1);
int3 = A3*cos(2*pi*fc3*t1);
int4(1:length(t1)) = 0;
random1 = randi([1 length(t1)]);random2 = randi([1 length(t1)]);random3 = randi([1 length(t1)]);
int4(random1) = 4;
int4(random2) = 2;
int4(random3) = 10;
interferenceSignals = int1 + int2 +int3 + int4;
modulatedSignalASK_noise_interf = modulatedSignalASKnoise + interferenceSignals;
figure(12);
plot(t1,modulatedSignalASK_noise_interf);
hold on;
plot(t1,int1,'b');
plot(t1,int2,'r');
plot(t1,int3,'m');
plot(t1,int4,'g');
title("ASK modulated OVER time, noise, interferences");

hold off;
figure(13);
plot(t1,modulatedSignalASK_noise_interf);
xlabel("Time [s]");
title("ASK modulated bits OVER TIME with noise and interference");

%%%ADDING ATTENUATION (DUE TO MEDIUM AND DISTANCE)
%%EACH MEDIUM HAS AN ATTENUATION CONSTANT alpha
%  exp(+-gamma*z)cos(i2pit)  =
%  exp(+-alpha*z)exp(+-iBetaz)cos(i2pit)
%%%%%%%%%%%%%%%%%%
taten = t1; %s
distMax = 10;
%distMax -> length(taten)
%zStep -> 1
zStep = 1*distMax/length(taten);
z = 0:zStep:distMax; % m
z = z(1:end-1);
alpha=0.3; %% 

c=3e8;
lambda = c/fc;
Beta = 2*pi/lambda; %aire
signalAten = exp(-alpha.*z).*exp(-i*Beta.*z).*modulatedSignalASK_noise_interf;
%signalAten = exp(-alpha.*z).*exp(-i*Beta.*z).*cos(2*pi*fcAt*taten);
figure(14);
plot(taten,signalAten);
xlabel("Time [s]");
title("ASK modulated bits OVER TIME with noise, interference and attenuation");
figure(15)
plot(z,signalAten);
xlabel("Distance [m]");
title("ASK modulated bits OVER DISTANCE with noise, interference and attenuation");

%%%ADDING PATH FADING (DUE TO OF THE MEDIUM TYPE AND THE EFFECTS ON SIGNAL SUCH AS REFLECTION, REFRACTION, SCATTER, ETC. AND DISTANCE)
%%WE CAN MODELED THIS BY PROPAGATION MODELS (i.g. Friis Model)%% Calculate
%%the path loss considering Free Space ans LoS considering a 1,2...,10m distance
%%and antennas gain of 3dBi, and fc=10kz., and Pt=20w o 43dBm;
%Lp = ...  


%%Recovering the signal (coherent demodulation)
%Multiplying
demodulatedSignalASK = signalAten.*s1;
demodulatedSignalASKfrec = abs(fftshift(fft(demodulatedSignalASK))); %DEP (PSD)
figure(16);
plot(freq,demodulatedSignalASKfrec,'m');
xlabel('Freq [Hz]');
title("Frequency Spectrum of ASK Demodulation with f_c= " + fc + "[Hz]");
%Filtering

demodulatedSignalASKfilter = filter(Hd,demodulatedSignalASK);
%demodulatedSignalASKfilter = demodulatedSignalASKfilter(1:end-16);
figure(17);
plot(tpulso,demodulatedSignalASKfilter);
xlabel('Time [s]');
title("ASK Demodulated bits in time");

demodulatedSignalASKfilterFrec = abs(fftshift(fft(demodulatedSignalASKfilter))); %DEP (PSD)
figure(18);
plot(freq,demodulatedSignalASKfilterFrec,'g');
xlabel('Freq [Hz]');
title("Frequency Spectrum of ASK bit DEmodulation Filtered");


%%Decision Making
%find the index of the timeClock Sample(half of period of pulse)
indexBit1 = find(tpulso==0.5);
indexBit2 = find(tpulso==1.5);
%%clock is set into the middle of each pulse (0.5s, 1.5s).

if(demodulatedSignalASKfilter(indexBit1) > 2000)
    bit(1) = 1;
else
    bit(1) = 0;
end

if(demodulatedSignalASKfilter(indexBit2) > 2000)
    bit(2) = 1;
else
    bit(2) = 0;
end

bit

% 1. Throughput probability = No. Frames Suc / ObsTime (No. Ranuras totales*TiempoCadaRanura)
% 
% 2. Error probability (Pe)  or BER = 1-No. FramesSuc / No. FramesTot
% 
% 3. Outage probability = (Nsuc – Ntot )/ Ntot
% 
% 4. Age of Information (AoI) = CMAT – PMAT (Current Message Arrival Time) – (Previous Message Arrival Time)
% 
% 5. Path Loss = Lp -> LoS es el modelo Friis. 
