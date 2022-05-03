close all
clear all
clc

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

%sendText = fopen('sendText.txt','w');
% for e=1:4557   
%     fprintf(sendText,'%c',H(e));
% end
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
T = length(H)-0.5;
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
N = fm*tmax;
freq = (-fm/2):(fm/N):(fm/2);
PS=fftshift(fft(pulseSignal));
figure(2)
plot(freq,abs(PS))
xlabel('Frecuencia [Hz]')
ylabel('Magnitud')
grid
title('Espectro de la señal de mensaje PS(w)')

N=100000;

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


%RECEPTOR 

%Demodulación coherente
demodulatedSignalASK = modulatedSignalASK.*c1;
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