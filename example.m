clear all;
clc;

%Supongamos una señal sinusoidal (campo electromagnético) propagándose por el aire a una
%frecuencia fc = 10Hz durante 5s;

fm = 44100; % frecuencia de muestreo [Hertz = muestras por seg];
Tm = 1/fm;


fc = 10000; %Hz
tmax=3; %tiempo en s
A = 20 ; %potencia de una torre celular, 43dBm = 20W 
t1 = 0: Tm :tmax;
c1 = A.*sin(2*pi*fc*t1);

figure(1);
plot(t1,c1);
xlabel('Time [s]');
title("Sine signal with f_c= " + fc + " [Hz]");

%%%Frequency analysis
figure(2);
fmax = fm/2;
N = fm*tmax; %muestras / s /s = muestras
f1 = -fmax : fm/N : fmax ;  %fs(muestras / s / muestras = 1 muestra /s)
C1 = abs(fftshift(fft(c1))); %DEP (PSD)
plot(f1,C1);
xlabel('Freq [Hz]');
title("Frequency Spectrum of sine signal with f_c= " + fc + " [Hz]");

%%Tarea: HAcer el ejercicio con diferentes frec y plot y saquen
%%conclusiones
%INDIVIDUAL
%PARA EL VIERNES 

%%%MODULACIÓN DIGITAL
pulsoBitUno(1:(length(t1)/2)+1)=A; %bit 1;
tpulso = t1;
pulsoBitCero(1:length(t1)/2)= 0; %Bit 0;


pulsos = horzcat(pulsoBitUno,pulsoBitCero);
figure(3);
plot(tpulso,pulsos,'r');
xlabel('Time [s]');
title("Bit Pulses 1 and 0");

%%
modulatedSignalASK = c1 .* pulsos;


figure(4);
plot(tpulso,modulatedSignalASK);
xlabel('Time [s]');
title("ASK bit modulation with f_c= " + fc + "[Hz]");


N = fm*tmax;
PULSOSfrec = abs(fftshift(fft(pulsos))); %DEP (PSD)
freq = (-fm/2):(fm/N):(fm/2);
figure(5);
plot(freq,PULSOSfrec,'g');
xlabel('Freq [Hz]');
title("Frequency Spectrum of 1 and 0 bits");

modulatedSignalASKfrec = abs(fftshift(fft(modulatedSignalASK))); %DEP (PSD)
freq = (-fm/2):(fm/N):(fm/2);
figure(6);
plot(freq,modulatedSignalASKfrec,'m');
xlabel('Freq [Hz]');
title("Frequency Spectrum of ASK bit modulation with f_c= " + fc + "[Hz]");


%%DEMODULEMOS LA SEÑAL POR MÉTODO COHERENTE

demodulatedSignalASK = modulatedSignalASK.*c1;
demodulatedSignalASKfrec = abs(fftshift(fft(demodulatedSignalASK))); %DEP (PSD)
figure(7);
plot(freq,demodulatedSignalASKfrec,'m');
xlabel('Freq [Hz]');
title("Frequency Spectrum of ASK Demodulation with f_c= " + fc + "[Hz]");


%Filtering
load filterLPfc10000.mat;
demodulatedSignalASKfilter = filter(Hd,demodulatedSignalASK);
%demodulatedSignalASKfilter = demodulatedSignalASKfilter(1:end-16);
figure(8);
plot(tpulso,demodulatedSignalASKfilter);
xlabel('Time [s]');
title("ASK Demodulated bits in time");

demodulatedSignalASKfilterFrec = abs(fftshift(fft(demodulatedSignalASKfilter))); %DEP (PSD)
figure(9);
plot(freq,demodulatedSignalASKfilterFrec,'g');
xlabel('Freq [Hz]');
title("Frequency Spectrum of ASK bit DEmodulation Filtered");


%%Decision Making
%find the index of the timeClock Sample(half of period of pulse)
indexBit1 = find(tpulso==0.75);
indexBit2 = find(tpulso==2.25);
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