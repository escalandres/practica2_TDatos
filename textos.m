close all
clear all
load FiltroLP400.mat
%%187 PALABRAS
%UTF-8 a Binario
fid = fopen('lyrics.txt');
b = fread(fid,'*uint8')';
fclose(fid);	

str = dec2bin(b,8);
disp(str);
fs = 44100; % frecuencia de muestreo [Hertz = muestras por seg];
Ts = 1/fs;
fc = 5010; %Hz
tmax=3; %tiempo en s
A = 20 ; %bit 1 en A=20W
t1 = 0: Ts :tmax;

%Creación del pulso unitario.
pulsoBitUno(1:(length(t1)/2)+1)=A; %bit 1;
tpulso = t1;
pulsoBitCero(1:length(t1)/2)= 0; %Bit 0;

pulseSignal=0;

for j=1:numel(str)/8
    for i=1:8
        if str(j,i)=='1'
            pulseSignal=horzcat(pulseSignal,pulsoBitUno);
        else
            pulseSignal=horzcat(pulseSignal,pulsoBitCero);
        end
    end
end
xmax=numel(str)*.01;
t2 = linspace(0,xmax,length(pulseSignal));

figure(1);
plot(t2,pulseSignal,'r');
xlabel('Time [s]');
axis([-0.1 xmax -0.1 5.1]);
title("Bit Pulses 1 and 0 with Bit 1 = 20, and Bit 0 = 0");

%Análisis en frecuencia de la señal moduladora (tren de pulsos)
N=100000;
w=linspace(-fs/2,fs/2,N)*2*pi;
PS=fftshift(fft(pulseSignal,N))*Ts;
figure(2)
plot(w/(2*pi),abs(PS))
xlabel('Frecuencia [Hz]')
ylabel('Magnitud')
grid
title('Espectro de la señal de mensaje PS(w)')

%Señal portadora
figure(3)
c1 = A*sin(2*pi*fc*t2);
subplot(211)
plot(t2,c1)
grid
title('Señal portadora c(t) with fc = 5010Hz and A = 20')

%Transformada de Fourier de la señal portadora
C1=fftshift(fft(c1,N))*Ts;
figure(3)
subplot(212)
plot(w/(2*pi),abs(C1))
grid
xlabel('Frecuencia [Hz]') ;
ylabel('Magnitud')
% axis([-400 400 0 0.09])
title('Espectro de la señal portadora C(w)')

%Señal modulada
modulatedSignalASK = c1 .* pulseSignal;
figure(4);
subplot(211)
plot(t2,modulatedSignalASK);
xlabel('Time [s]');
title("ASK bit modulation with f_c= " + fc + "[Hz]");

%Transformada de Fourier de la señal Modulada
modulatedSignalASKfrec=fftshift(fft(modulatedSignalASK,N))*Ts;
figure(4)
subplot(212)
plot(w/(2*pi),abs(modulatedSignalASKfrec))
xlabel('Freq [Hz]');
ylabel('Magnitud')
grid
title("Frequency Spectrum of ASK bit modulation with f_c= " + fc + "[Hz]");


%RECEPTOR 

%Demodulación coherente
demodulatedSignalASK = modulatedSignalASK.*c1;
demodulatedSignalASKfrec = abs(fftshift(fft(demodulatedSignalASK,N))); %DEP (PSD)
figure(5);
plot(w,demodulatedSignalASKfrec,'m');
xlabel('Freq [Hz]');
title("Frequency Spectrum of ASK Demodulation with f_c= " + fc + "[Hz]");

%Filtro
demodulatedSignalASKfilter = filter(Hd,demodulatedSignalASK);
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
