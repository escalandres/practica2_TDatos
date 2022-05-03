close all
clear all
clc
%Practica 2.2: Ruido
%Parte 1: conversion del texto a binario

%126 palabras
text = 'A long, long time ago. I can still remember how that music used to make me smile. And I knew if I had my chance that I could make those people dance. And maybe they would be happy for a while. But February made me shiver. With every paper I would deliver. Bad news on the doorstep. I could not take one more step. I cannot remember if I cried. When I read about his widowed bride. But something touched me deep inside. The day the music died. So bye-bye, Miss American Pie. Drove my Chevy to the levee, but the levee was dry. And them good old boys were drinking whiskey and rye. Singin This will be the day that I die. This will be the day that I die';
%text = 'A long, long time ago.';
binary = dec2bin(text);
%str = bin2dec(binary);

%Parte 2: convertir binary a un vector (1x4557) 1 row and 4557 colums

%651*7 = 4557 --651 letras de 7 bits cada una
limit = length(text)*7;
%Vector de 1x4557, donde se guardaran todos los datos binarios de binary
H = zeros ([1 limit]);
count = 1;
for i = 1:length(text)
    for j = 1:7
%         h(count) = str2num(binary(i,j));
        H(count) = str2double(binary(i,j));
        count=count+1;
    end
end

sendText = fopen('sendText.txt','w');
for e=1:4557   
    fprintf(sendText,'%c',H(e));
end
disp('Binary convertion completed')

disp('binary sended')
disp(H)

A = 20; %W potencia de una torre celular, 43dBm = 20W 
fm = 44100; % frecuencia de muestreo [Hertz = muestras por seg];
Tm = 1/fm;
fc = 5010; %Hz
tmax=3; %tiempo en s
t1 = 0: Tm :tmax;
disp('ASK modulation starting')
%%%MODULACIÓN DIGITAL
T = length(H);
pulsoBitUno(1:(length(t1)/T)+1)=A; %bit 1;
tpulso = t1;
pulsoBitCero(1:length(t1)/T)=-A; %Bit 0;

pulseSignal=0;

for k = 1:length(H)
    if H(k)==1
        pulseSignal=horzcat(pulseSignal,pulsoBitUno);
    else
        pulseSignal=horzcat(pulseSignal,pulsoBitCero);
    end
end

%xmax=numel(binary)*.01;
%t2 = linspace(0,xmax,length(pulseSignal));
add = 3.0000;

for k = length(t1):length(t1)+(length(pulseSignal)-length(t1))
    t1(k) = add; 
end

figure(1);
plot(t1,pulseSignal,'r');
xlabel('Time [s]');
axis([-0.1 3.1 -A-0.2 A+0.2]);
title("Bit Pulses 1 and 0 with Bit 1 = "+A+", and Bit 0 = "+(-A));
disp('ASK modulation completed')

%Análisis en frecuencia de la señal moduladora (tren de pulsos)
N=100000;
w=linspace(-fm/2,fm/2,N)*2*pi;
PS=fftshift(fft(pulseSignal,N))*Tm;
figure(2)
plot(w/(2*pi),abs(PS))
xlabel('Frecuencia [Hz]')
ylabel('Magnitud')
grid
title('Espectro de la señal de mensaje PS(w)')

%Señal portadora
figure(3)
c1 = A*sin(2*pi*fc*t1);
subplot(211)
plot(t1,c1)
grid
title("Señal portadora c(t) with fc = "+fc+" and A = "+A)

%Transformada de Fourier de la señal portadora
C1=fftshift(fft(c1,N))*Tm;
figure(3)
subplot(212)
plot(w/(2*pi),abs(C1))
grid
xlabel('Frecuencia [Hz]') ;
ylabel('Magnitud')
% axis([-400 400 0 0.09])
title('Espectro de la señal portadora C(w)')

disp('Carry signal Fourier transform completed')

%Señal modulada
modulatedSignalASK = c1 .* pulseSignal;
figure(4);
subplot(211)
plot(t1,modulatedSignalASK);
xlabel('Time [s]');
title("ASK bit modulation with f_c= " + fc + "[Hz]");

%Transformada de Fourier de la señal Modulada
modulatedSignalASKfrec=fftshift(fft(modulatedSignalASK,N))*Tm;
figure(4)
subplot(212)
plot(w/(2*pi),abs(modulatedSignalASKfrec))
xlabel('Freq [Hz]');
ylabel('Magnitud')
grid
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













%RECEPTOR 

%Demodulación coherente
demodulatedSignalASK = signalAten.*c1;
demodulatedSignalASKfrec = abs(fftshift(fft(demodulatedSignalASK,N))); %DEP (PSD)
figure(5);
plot(w,demodulatedSignalASKfrec,'m');
xlabel('Freq [Hz]');
title("Frequency Spectrum of ASK Demodulation with f_c= " + fc + "[Hz]");

%Filtro
%load filterLPK5010.mat % 
%demodulatedSignalASKfilter = filter(Hd,demodulatedSignalASK);
load filterLP126W.mat % 
demodulatedSignalASKfilter = filter(Fd,demodulatedSignalASK);
%demodulatedSignalASKfilter = demodulatedSignalASKfilter(1:end-16);
figure(6);
plot(demodulatedSignalASKfilter);
xlabel('Time [s]');
title("ASK Demodulated bits in time");

demodulatedSignalASKfilterFrec = abs(fftshift(fft(demodulatedSignalASKfilter))); %DEP (PSD)
figure(7);
plot(demodulatedSignalASKfilterFrec,'g');
xlabel('Freq [Hz]');
title("Frequency Spectrum of ASK bit DEmodulation Filtered");

%%Decision Making
%find the index of the timeClock Sample(half of period of pulse)
% indexBit1 = find(tpulso==0.75);
% indexBit2 = find(tpulso==2.25);
%%clock is set into the middle of each pulse (0.5s, 1.5s).

%Caso 22 palabras
TRecuperado = length(demodulatedSignalASKfilter)/length(H);
limitRecuperado = length(demodulatedSignalASKfilter)/TRecuperado;
limitRecuperado1 = round(limitRecuperado);
bitsRecuperados = zeros([1 limitRecuperado1]);

kBitRecuperados = 1;
for k = 570:round(TRecuperado):length(demodulatedSignalASKfilter)
    if(demodulatedSignalASKfilter(k) > 2000)
        bitsRecuperados(kBitRecuperados) = 1;
    else
        bitsRecuperados(kBitRecuperados) = 0;
    end
    kBitRecuperados = kBitRecuperados + 1;
end
disp('Bits recuperados')
%disp(bitsRecuperados)

bitsRecovered = zeros([1 length(H)]);
Rcount = 1;
for k=1:4559
    bitsRecovered(Rcount) = bitsRecuperados(k);
    Rcount = Rcount+1;
end

bitsRecovered(32) = [];
bitsRecovered(82) = [];
disp(bitsRecovered)

%Conversion del binario recuperado en texto
filas = round(length(bitsRecovered)/7);
str = zeros([filas 7]);
bitCount = 1;
for i = 1:filas
    for j = 1:7
        str(i,j) = bitsRecovered(bitCount);
        bitCount=bitCount+1;
    end
end
pA = string(str);
pru = zeros([filas 1]);

for k=1:filas
    pru(k) = bin2dec(pA(k,1)+pA(k,2)+pA(k,3)+pA(k,4)+pA(k,5)+pA(k,6)+pA(k,7));
end

strArray = pru;
textRecuperado = native2unicode(strArray);
disp('Texto recuperado:')
disp(textRecuperado.')