close all
clear all
clc

%Parte 1: conversion del texto a binario

%126 palabras
%text = 'A long, long time ago. I can still remember how that music used to make me smile. And I knew if I had my chance that I could make those people dance. And maybe they would be happy for a while. But February made me shiver. With every paper I would deliver. Bad news on the doorstep. I could not take one more step. I cannot remember if I cried. When I read about his widowed bride. But something touched me deep inside. The day the music died. So bye-bye, Miss American Pie. Drove my Chevy to the levee, but the levee was dry. And them good old boys were drinking whiskey and rye. Singin This will be the day that I die. This will be the day that I die';
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
%H = FCRC(H);

% bS = "";
% for i=1:length(H)
%     a = string(H(i));
%     %bS = horzcat(bS,a);
%     bS = bS+a;
% end
% 
% smessage="0b"+bS+"u32";
% %message = 0b1101100111011010u32;
% message = str2num(smessage);
% messageLength = length(H);
% divisor = 0b1101u32;
% divisorDegree= 3;
% 
% divisor = bitshift(divisor,messageLength-divisorDegree-1);
% dec2bin(divisor)
% 
% divisor = bitshift(divisor,divisorDegree);
% remainder = bitshift(message,divisorDegree);
% dec2bin(divisor)
% 
% dec2bin(remainder)
% 
% for k=1:messageLength
% 	if bitget(remainder,messageLength+divisorDegree)
% 		remainder = bitxor(remainder,divisor);
% 	end
% 	remainder = bitshift(remainder,1);
% end
% 
% CRC_check_value = bitshift(remainder,-messageLength);
% FCS = dec2bin(CRC_check_value);
% H(length(H)+1)=FCS(1);
% H(length(H)+1)=FCS(2);
% H(length(H)+1)=FCS(3);


disp('Binary convertion completed')

disp('binary sended')
% disp(H)

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
tpulso = t1;

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

%%%ADDING NOISE AWGN
SNR=100; %%dBs
modulatedSignalASKnoise = awgn(modulatedSignalASK,SNR,'measured');
figure(6);
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
figure(7);
plot(t1,modulatedSignalASK_noise_interf);
hold on;
plot(t1,int1,'b');
plot(t1,int2,'r');
plot(t1,int3,'m');
plot(t1,int4,'g');
title("ASK modulated OVER time, noise, interferences");

hold off;
figure(8);
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
figure(9);
plot(taten,signalAten);
xlabel("Time [s]");
title("ASK modulated bits OVER TIME with noise, interference and attenuation");
figure(10)
plot(z,signalAten);
xlabel("Distance [m]");
title("ASK modulated bits OVER DISTANCE with noise, interference and attenuation");

%%%ADDING PATH FADING (DUE TO OF THE MEDIUM TYPE AND THE EFFECTS ON SIGNAL SUCH AS REFLECTION, REFRACTION, SCATTER, ETC. AND DISTANCE)
%%WE CAN MODELED THIS BY PROPAGATION MODELS (i.g. Friis Model)%% Calculate
%%the path loss considering Free Space ans LoS considering a 1,2...,10m distance
%%and antennas gain of 3dBi, and fc=10kz., and Pt=20w o 43dBm;
%Lp = ...  
% figure(11)
% Ptx = 43; %dBm, Potencia de la torre celular
% Gtx = 3; %dBi, Ganancia de la torre celular
% Grx = 3; %dBi, Ganancia del receptor
% frec = 10e3; %MHZ, Frecuencia de la Torre celular
% lamba = c/frec;
% R = 1:10; %m, Distancia entre el receptor y la Torre celular
% Prx = Ptx+Gtx+Grx+20*log(lamba)-20*log(4*pi*R);
% Plost = Ptx - Prx;
% plot(R,Plost,"LineWidth",2,"Color",'b')
% xlabel('Distancia R[m]')
% ylabel('Perdida[dBm]')
% title('Perdida por trayectoria Friis')



%RECEPTOR 

%%Recovering the signal (coherent demodulation)
%Multiplying
demodulatedSignalASK = signalAten.*c1;
figure(11)
plot(tpulso,demodulatedSignalASK)
demodulatedSignalASKfrec = abs(fftshift(fft(demodulatedSignalASK))); %DEP (PSD)
figure(12);
plot(freq,demodulatedSignalASKfrec,'m');
xlabel('Freq [Hz]');
title("Frequency Spectrum of ASK Demodulation with f_c= " + fc + "[Hz]");




%Filtro
load filterLP126Noise.mat % 
%demodulatedSignalASKfilter = filter(Hd,demodulatedSignalASK);
%load filterLP126W.mat % 
demodulatedSignalASKfilter = filter(Fd,demodulatedSignalASK);
%demodulatedSignalASKfilter = demodulatedSignalASKfilter(1:end-16);
figure(13);
plot(demodulatedSignalASKfilter);
xlabel('Time [s]');
title("ASK Demodulated bits in time");

demodulatedSignalASKfilterFrec = abs(fftshift(fft(demodulatedSignalASKfilter))); %DEP (PSD)
figure(14);
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
    if(demodulatedSignalASKfilter(k) > 90)
        bitsRecuperados(kBitRecuperados) = 1;
    else
        bitsRecuperados(kBitRecuperados) = 0;
    end
    kBitRecuperados = kBitRecuperados + 1;
end
disp('Bits recuperados')
%disp(bitsRecuperados)
bitsRecovered=bitsRecuperados;

% smessage="0b"+bitsRecovered+"u32";
% message = 0b1101100111011010u32;
% message = str2num(smessage);
% messageLength = length(bitStream);
% divisor = 0b1101u32;
% divisorDegree= 3;
% 
% divisor = bitshift(divisor,messageLength-divisorDegree-1);
% dec2bin(divisor)
% 
% divisor = bitshift(divisor,divisorDegree);
% remainder = bitshift(message,divisorDegree);
% dec2bin(divisor)
% 
% dec2bin(remainder)
% 
% for k=1:messageLength
% 	if bitget(remainder,messageLength+divisorDegree)
% 		remainder = bitxor(remainder,divisor);
% 	end
% 	remainder = bitshift(remainder,1);
% end
% 
% CRC_check_value = bitshift(remainder,-messageLength);
% FCS = dec2bin(CRC_check_value);

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
textRecuperado=textRecuperado.';
disp(textRecuperado)

%1. PATH LOSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ptx = 43;
Gtx = 12;
Grx = 3;
c = 3e8;
distMax = 10;
d = 0:distMax;

P = Ptx + Gtx + Grx + 20*log(c/fc) - 20*log(4*pi*d);
figure(17)
plot(d,P);
xlabel("Distance [m]");
ylabel('Potencia');
grid
title("ASK modulated bits OVER DISTANCE with noise, interference and attenuation");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mensaje = strArray.';
mensaje = char(mensaje);
mensaje1 = text;
mensaje1 = char(mensaje1);


%2. THROUGHPUT
framesFallidos = 0;

framesTotales = length(mensaje);
for i = 1:length(mensaje)
    if strcmp(mensaje(i), mensaje1(i)) == 0
        framesFallidos = framesFallidos + 1;
    end
end
t2 = tpulso;
framesExitosos = framesTotales - framesFallidos;
tiempo = (length(t2)-1) / 11000;
throuhput = framesExitosos/tiempo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3. ERROR POSIBILITY
FER = framesFallidos/framesTotales;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4. OUTTAGE POSIBILITY
OUTP = (framesExitosos-framesTotales)/framesTotales;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%5. AGE OF INFORMATION
AOI = ((tiempo/framesTotales)*2 - tiempo/framesTotales);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CRC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aux = str2num(str(1,:))

% for i = 1:length(str)
%     listaSE(i) = bin2dec(str(i,:))
% end
%     listaCE = mensaje 
listaCE = strArray.';
% for k=1:length(binary)
%     listaSE(k) = bin2dec(binary);
% end
listaSE = strArray.';
errores = 0;
for i = 1:length(binary) 
    message = listaSE(i);
    messageLength = 8;
    divisor = 0b1101u32;
    divisorDegree = 3;
    
    divisor = bitshift(divisor,messageLength-divisorDegree-1);
    dec2bin(divisor)
    
    divisor = bitshift(divisor,divisorDegree);
    remainder = bitshift(message,divisorDegree);
    dec2bin(divisor)
    
    dec2bin(remainder)
    
    for k = 1:messageLength
        if bitget(remainder,messageLength+divisorDegree)
            remainder = bitxor(remainder,divisor);
        end
        remainder = bitshift(remainder,1);
    end
    
    CRC_check_value = bitshift(remainder,-messageLength);
    dec2bin(CRC_check_value)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %COMPROBACION DEL MENSAJE
    
    % remainder = bitset(remainder,[2 3]);
    % remainder = bitset(remainder,6);
    
    for k = 1:messageLength
        remainder = bitshift(listaCE(i),divisorDegree);
        remainder = bitor(remainder,CRC_check_value);
        if bitget(remainder,messageLength+divisorDegree)
            remainder = bitxor(remainder,divisor);
        end
        remainder = bitshift(remainder,1);
    end
    if remainder == 0
        disp('Message is error free.')
    
    else
        disp('Message contains errors.')
        errores = errores + 1;
    end
end 