close all
clear all
clc

%Parte 1: conversion del texto a binario

%126 palabras
%text = 'A long, long time ago. I can still remember how that music used to make me smile. \nAnd I knew if I had my chance that I could make those people dance. \nAnd maybe they would be happy for a while. But February made me shiver. \nWith every paper I would deliver. Bad news on the doorstep. I could not take one more step. \nI cannot remember if I cried. When I read about his widowed bride. But something touched me deep inside. \nThe day the music died. So bye-bye, Miss American Pie. Drove my Chevy to the levee, but the levee was dry. \nAnd them good old boys were drinking whiskey and rye. Singin This will be the day that I die. This will be the day that I die';
%text = 'A long, long time ago.';
fid = fopen('AmericanPieLyrics.txt');
b = fread(fid,'*uint8')';
fclose(fid);
text = b;
binary = dec2bin(text,8);
%str = bin2dec(binary);

%Parte 2: convertir binary a un vector (1x4557) 1 row and 4557 colums

%651*7 = 4557 --651 letras de 7 bits cada una
limit = length(text)*8;
%Vector de 1x4557, donde se guardaran todos los datos binarios de binary
H = zeros ([1 limit]);
count = 1;
for i = 1:length(text)
    for j = 1:8
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

%%%MODULACIÓN DIGITAL
Fbit = 100;
Tbit = 1/Fbit;
disp('Polar NRZ codification started')
bitStream=TestPNRZ(H,Tbit);
bitStream(1)=[];
bitStream(length(bitStream))=[];
A = 20; %W potencia de una torre celular, 43dBm = 20W 
for k=1:length(bitStream)
    if bitStream(k)==1
        bitStream(k)=A;
    else
        bitStream(k)=-A;
    end
end
disp('Polar NRZ codification finished')

fm = 44100; % frecuencia de muestreo [HertbitStream = muestras por seg];
Tm = 1/fm;
fc = 5010; %Hz
tmax=(3*length(bitStream))/132301; %tiempo en s
t1 = 0: Tm :tmax;
t1(length(t1)+1)=10.4366;
t1(length(t1)+1)=10.4366;
t1(length(t1)+1)=10.4366;

figure(2);
plot(t1,bitStream,'r');
xlabel('Time [s]');
axis([-0.1 tmax+1 -A-2 A+2]);
title("Bit Pulses 1 and 0 with Bit 1 = "+A+", and Bit 0 = "+(-A));
disp('ASK modulation completed')

%Señal portadora
figure(3)
c1 = A*sin(2*pi*fc*t1);
subplot(211)
plot(t1,c1)
grid
title("Señal portadora c(t) with fc = "+fc+" and A = "+A)

%Transformada de Fourier de la señal portadora
figure(3)
fmax = fm/2;
N = fm*tmax; %muestras / s /s = muestras
f1 = -fmax : fm/N : fmax ;  %fs(muestras / s / muestras = 1 muestra /s)
f1(length(f1)+1)=2.2050;
f1(length(f1)+1)=2.2050;
f1(length(f1)+1)=2.2050;
C1 = abs(fftshift(fft(c1))); %DEP (PSD)
subplot(212)
plot(f1,C1);
xlabel('Frecuencia [Hz]') ;
ylabel('Magnitud')
title("Frequency Spectrum of sine signal with f_c= " + fc + " [Hz]");


% axis([-400 400 0 0.09])
title("Espectro de la señal portadora C(w) con f_c= " + fc + " [Hz]")
disp('Carry signal Fourier transform completed')

%Señal modulada
modulatedSignalASK = c1 .* bitStream;
figure(4);
subplot(211)
plot(t1,modulatedSignalASK);
xlabel('Time [s]');
axis([-1 tmax+1 -430 430])
title("ASK bit modulation with f_c= " + fc + "[Hz]");

%Transformada de Fourier de la señal Modulada
modulatedSignalASKfrec = abs(fftshift(fft(modulatedSignalASK))); %DEP (PSD)
freq = (-fm/2):(fm/N):(fm/2);
freq(length(freq)+1)=2.2050;
freq(length(freq)+1)=2.2050;
freq(length(freq)+1)=2.2050;
figure(4);
subplot(212)
plot(freq,modulatedSignalASKfrec,'m');
xlabel('Freq [Hz]');
ylabel('Magnitud')
title("Frequency Spectrum of ASK bit modulation with f_c= " + fc + "[Hz]");

N = fm*tmax;
PULSOSfrec = abs(fftshift(fft(bitStream))); %DEP (PSD)
% freq = (-fm/2):(fm/N):(fm/2);
figure(5);
plot(freq,PULSOSfrec,'g');
xlabel('Freq [Hz]');
title("Frequency Spectrum of 1 and 0 bits");

%RECEPTOR 

%Demodulación coherente
demodulatedSignalASK = modulatedSignalASK.*c1;
demodulatedSignalASKfrec = abs(fftshift(fft(demodulatedSignalASK))); %DEP (PSD)
figure(6);
plot(freq,demodulatedSignalASKfrec,'m');
xlabel('Freq [Hz]');
title("Frequency Spectrum of ASK Demodulation with f_c= " + fc + "[Hz]");

%Filtro
%load filterLPK5010.mat % 
%demodulatedSignalASKfilter = filter(Hd,demodulatedSignalASK);
load filterLP126W.mat % 
demodulatedSignalASKfilter = filter(Fd,demodulatedSignalASK);
%demodulatedSignalASKfilter = demodulatedSignalASKfilter(1:end-16);
figure(7);
plot(demodulatedSignalASKfilter);
xlabel('Time [s]');
title("ASK Demodulated bits in time");

demodulatedSignalASKfilterFrec = abs(fftshift(fft(demodulatedSignalASKfilter))); %DEP (PSD)
figure(8);
plot(demodulatedSignalASKfilterFrec,'g');
xlabel('Freq [Hz]');
title("Frequency Spectrum of ASK bit DEmodulation Filtered");

%%Decision Making

%Caso 22 palabras -------------------
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

% bitsRecovered = zeros([1 length(H)]);
% Rcount = 1;
% for k=1:4559
%     bitsRecovered(Rcount) = bitsRecuperados(k);
%     Rcount = Rcount+1;
% end

% bitsRecovered(32) = [];
% bitsRecovered(82) = [];
% disp(bitsRecovered)
bitsRecovered=bitsRecuperados;

%Conversion del binario recuperado en texto
filas = round(length(bitsRecovered)/8);
str = zeros([filas 8]);
bitCount = 1;
for i = 1:filas
    for j = 1:8
        str(i,j) = bitsRecovered(bitCount);
        bitCount=bitCount+1;
    end
end
pA = string(str);
pru = zeros([filas 1]);

for k=1:filas
    pru(k) = bin2dec(pA(k,1)+pA(k,2)+pA(k,3)+pA(k,4)+pA(k,5)+pA(k,6)+pA(k,7)+pA(k,8));
end

strArray = pru;
textRecuperado = native2unicode(strArray);
disp('Texto recuperado:')
disp(textRecuperado.')