
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DISEÑO DEL FILTRO KAISER TIPO VENTANA
%CON FRECUENCIA DE MUESTRO (fs=11000)
%FPASS = 50
%Y FRECUENCIA DE CORTE = 100
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hd = 0
Fs = 11000;  % Sampling Frequency
Fpass = 50;               % Passband Frequency
Fstop = 100;              % Stopband Frequency
Dpass = 0.057501127785;  % Passband Ripple
Dstop = 0.0001;          % Stopband Attenuation
flag  = 'scale';         % Sampling Flag
[N,Wn,BETA,TYPE] = kaiserord([Fpass Fstop]/(Fs/2), [1 0], [Dstop Dpass]);
% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
Hd = dfilt.dffir(b);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TRANSMISOR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('texto2.txt');
b = fread(fid);
fclose(fid);	
str = dec2bin(b,8);
disp(str);

fs = 11000;
Ts = 1/fs;
fc = 5000; 
A = 20; %20W transmitted power (43dBm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREACION DE LOS BITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pulsoBitUno(1:fs) = A; 
pulsoBitCero(1:fs) = 0; 
pulseSignal=0;

for j=1:numel(str)/8
    for i=1:8
        if str(j,i) == '1'
            pulseSignal=horzcat(pulseSignal,pulsoBitUno);
        else
            pulseSignal=horzcat(pulseSignal,pulsoBitCero);
        end
    end

xmax=numel(str)*1;
t2 = linspace(0,xmax,length(pulseSignal));

figure(1);
plot(t2,pulseSignal,'r');
xlabel('Time [s]');
ylabel('Magnitude')
axis([-0.1 xmax -0.1 20.1]);
grid
title("Bit Pulses 1 and 0");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALISIS EN FRECUENCIA DE LA SEÑAL MODULADORA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=100000;
w=linspace(-fs/2,fs/2,N)*2*pi;
PS=fftshift(fft(pulseSignal,N))*Ts;

figure(2)
plot(w/(2*pi),abs(PS))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
axis([-100 100 -.1 43])
grid
title('Espectro de la señal de mensaje PS(w)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SEÑAL PORTADORA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
c1 = A*sin(2*pi*t2*fc);
subplot(211)
plot(t2,c1)
xlabel('Frequency [Hz]')
ylabel('Magnitude')
axis([-0.1 xmax -20.1 20.1]);
grid
title('Señal portadora c(t)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRANSFORMADA DE LA PORTADORA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C1=fftshift(fft(c1,N))*Ts;
figure(3)
subplot(212)
plot(w/(2*pi),abs(C1))
xlabel('Frequency [Hz]') 
ylabel('Magnitude')
axis([-5500 5500 0 65])
grid
title('Espectro de la señal portadora C(w)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SEÑAL MODULADA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modulatedSignalASK = c1 .* pulseSignal;
figure(4);
subplot(211)
plot(t2,modulatedSignalASK);
xlabel('Time [s]');
ylabel('Magnitude')
axis([0 xmax -400 400])
grid
title("ASK bit modulation with f_c= " + fc + "[Hz]");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRANSFORMADA DE LA SEÑAL MODULADA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modulatedSignalASKfrec=fftshift(fft(modulatedSignalASK,N))*Ts;
figure(4)
subplot(212)
plot(w/(2*pi),abs(modulatedSignalASKfrec))
xlabel('Frequency [Hz]')
ylabel('Magnitude')
axis([-5500 5000 0 360])
grid
title("Frequency Spectrum of ASK bit modulation with f_c= " + fc + "[Hz]");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%RECEPTOR 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEMODULACION COHERENTE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
demodulatedSignalASK = modulatedSignalASK.*c1;
demodulatedSignalASKfrec = fftshift(fft(demodulatedSignalASK,N)); 
figure(5);
plot(w,abs(demodulatedSignalASKfrec),'m');
xlabel('Frequency [Hz]');
ylabel('Magnitude')
grid
title("Frequency Spectrum of ASK Demodulation with f_c= " + fc + "[Hz]");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FILTRO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
demodulatedSignalASKfilter = filter(Hd,demodulatedSignalASK);
figure(6);
plot(t2, demodulatedSignalASKfilter);
xlabel('Time [s]');
ylabel('Magnitude')
axis([0 xmax -400 4350])
grid
title("ASK Demodulated bits in time");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SEÑAL DEMODULADA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
demodulatedSignalASKfilterFrec = abs(fftshift(fft(demodulatedSignalASKfilter, N))); %DEP (PSD)
figure(7);
plot(w,demodulatedSignalASKfilterFrec,'g');
xlabel('Frequency [Hz]');
ylabel('Magnitude')
axis([-100 100 0 9e7])
grid
title("Frequency Spectrum of ASK bit Demodulation Filtered");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RECUPERACION DE LA SEÑAL DE MENSAJE ENVIADA EN BINARIO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mensajeRecuperado = 0
aux = 0
a = 1
b = 1
contador = 1
for j = 1:xmax
    if demodulatedSignalASKfilter((11000*j)-5500-1) > 30 && j == 1
        aux(a,b) = 1
    elseif demodulatedSignalASKfilter((11000*j)-5500-1) < 30  && j == 1
        aux(a,b) = 0
    elseif demodulatedSignalASKfilter((11000*j)-5500) > 30  
        aux(a,b) = 1
    elseif demodulatedSignalASKfilter((11000*j)-5500) < 30  
        aux(a,b) = 0
    end
    if j == 1
        contador = contador + 11000 -1
    else
        contador = contador + 11000 
    end
    if b == 8
        a = a + 1
        b = 1
    else
        b = b +1
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DECODIFICACIÓN DE LA SEÑAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aux2 = []
mensaje1 = []
contador = 1
for j = 1:numel(aux)/8
    for i = 1:8
        aux2(i) = aux(j,i)
    end
    mensaje1(contador) = bin2dec(num2str(aux2))
    contador = contador + 1
end
mensaje1 = char(mensaje1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SEGUNDA PARTE DE LA PRACTICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SEÑALES PORTADORAS QUE SERVIRAN DE INTERFERENCIA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = A/3*sin(2*pi*3*fc*t2);
s2 = A*sin(2*pi*5*fc*t2+pi);

figure(8)
subplot(211)
plot(t2, s1)
xlabel('Time [s]');
ylabel('Magnitude')
axis([-0.1 xmax -6.7 6.7]);
grid
title 'Interference Signal 1 (A/3*sin(2π*2f_c*t)'

subplot(212)
plot(t2, s2)
xlabel('Time [s]');
ylabel('Magnitude')
axis([-0.1 xmax -20.1 20.1]);
grid
title 'Interference Signal 2 (A*sin(2π*3f_c*t + π)'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OBTENIENDO LA SEÑAL PORTADORA A TRABAJAR CON INTERFERENCIA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sTotal = c1 + s1 + s2;
% 
% figure(9)
% subplot(211)
% plot(t2, sTotal)
% xlabel('Time [s]');
% ylabel('Magnitude')
% %axis([-0.1 xmax -6.7 6.7]);
% hold on
% plot(t2, s2, 'g')
% hold on
% plot(t2, s1, 'r')
% grid
% title 'Carrier signal with interference'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ESPECTRO DE LA SEÑAL PORTADORA CON INTERFERENCIA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STOTAL=fftshift(fft(sTotal,N))*Ts
% subplot(212)
% plot(w, abs(STOTAL))
% xlabel('Frequency [Hz]');
% ylabel('Magnitude')
% %axis([-0.1 xmax -20.1 20.1]);
% grid
% title 'Carrier signal spectrum with interference'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AGREGANDO RUIDO AWGN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR = 30; %%dBs
modulatedSignalASKnoise = awgn(modulatedSignalASK,SNR,'measured');
figure(9);
plot(t2,modulatedSignalASKnoise);
xlabel('Time [s]');
ylabel('Magnitude');
grid
title("ASK signal with awgn noise: " + SNR +"dB");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AGREGANDO RUIDO AWGN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A1 = 20; 
A2 = 20; 
A3 = 20;
fc1 = 10000; 
fc2 = 900; 
fc3 = 500;
int1 = A1*cos(2*pi*fc1*t2);
int2 = A2*cos(2*pi*fc2*t2);
int3 = A3*cos(2*pi*fc3*t2);
int4(1:length(t2)) = 0;
random1 = randi([1 length(t2)]);
random2 = randi([1 length(t2)]);
random3 = randi([1 length(t2)]);
int4(random1) = 5;
int4(random2) = 7;
int4(random3) = 12;
interferenceSignals = int1 + int2 +int3 + int4;
modulatedSignalASK_noise_interf = modulatedSignalASKnoise + interferenceSignals;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SEÑAL MODULADA CON RUIDO E INTERFERENCIA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(10);
plot(t2,modulatedSignalASK_noise_interf);
hold on;
plot(t2,int1,'b');
plot(t2,int2,'r');
plot(t2,int3,'m');
plot(t2,int4,'g');
grid
xlabel('Time [s]');
ylabel('Magnitude');
grid
title("ASK modulated OVER time, noise, interferences");

figure(11);
plot(t2,modulatedSignalASK_noise_interf);
xlabel("Time [s]");
ylabel('Magnitude');
grid
title("ASK modulated bits OVER TIME with noise and interference");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ATENUACION DE LA SEÑAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
taten = t2; 
distMax = 100;
zStep = 1*distMax/length(taten);
z = 0:zStep:distMax; 
z = z(1:end-1);
alpha=0.01; 
c=3e8;
lambda = c/fc;
Beta = 2*pi/lambda; 
signalAten = exp(-alpha.*z).*exp(-i*Beta.*z).*modulatedSignalASK_noise_interf;

figure(12);
plot(taten,signalAten);
xlabel("Time [s]");
ylabel('Magnitude');
grid
title("ASK modulated bits OVER TIME with noise, interference and attenuation");

figure(13)
plot(z,signalAten);
xlabel("Distance [m]");
ylabel('Magnitude');
grid
title("ASK modulated bits OVER DISTANCE with noise, interference and attenuation");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COHERENT DEMODULATION WITH NOISE AND INTERFERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
demodulatedSignalASK = signalAten.*c1;
demodulatedSignalASKfrec = fftshift(fft(demodulatedSignalASK, N)) 
figure(14);
plot(w/2*pi,abs(demodulatedSignalASKfrec),'m');
xlabel('Frequency [Hz]');
ylabel('Magnitude');
grid
title("Frequency Spectrum of ASK Demodulation with f_c= " + fc + "[Hz]");

%FILTRO
demodulatedSignalASKfilter = filter(Hd,demodulatedSignalASK);
figure(15);
plot(t2,demodulatedSignalASKfilter);
xlabel('Time [s]');
ylabel('Magnitude');
grid
title("ASK Demodulated bits in time");

%ESPECTRO DE LA SEÑAL OBTENIDA
demodulatedSignalASKfilterFrec = abs(fftshift(fft(demodulatedSignalASKfilter, N)));
figure(16);
plot(w/2*pi,demodulatedSignalASKfilterFrec,'g');
xlabel('Frequency [Hz]');
ylabel('Magnitude');
grid
title("Frequency Spectrum of ASK bit Demodulation Filtered");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RECUPERACION DE LA SEÑAL DE MENSAJE ENVIADA EN BINARIO (CON RUIDO E
%INTERFERENCIA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mensajeRecuperado = 0
aux = 0
a = 1
b = 1
contador = 1
for j = 1:xmax
    if demodulatedSignalASKfilter((11000*j)-5500-1) > 1500 && j == 1
        aux(a,b) = 1
    elseif demodulatedSignalASKfilter((11000*j)-5500-1) < 1500  && j == 1
        aux(a,b) = 0
    elseif demodulatedSignalASKfilter((11000*j)-5500) > 1500  
        aux(a,b) = 1
    elseif demodulatedSignalASKfilter((11000*j)-5500) < 1500  
        aux(a,b) = 0
    end
    if j == 1
        contador = contador + 11000 -1
    else
        contador = contador + 11000 
    end
    if b == 8
        a = a + 1
        b = 1
    else
        b = b +1
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DECODIFICACIÓN DE LA SEÑAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aux2 = []
mensaje = []
contador = 1
for j = 1:numel(aux)/8
    for i = 1:8
        aux2(i) = aux(j,i)
    end
    mensaje(contador) = bin2dec(num2str(aux2))
    contador = contador + 1
end
mensajeRecuperado = char(mensaje)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. PATH LOSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ptx = 43
Gtx = 12
Grx = 3
c = 3e8
d = 0:distMax

P = Ptx + Gtx + Grx + 20*log(c/fc) - 20*log(4*pi*d)

figure(17)
plot(d,P);
xlabel("Distance [m]");
ylabel('Potencia');
grid
title("ASK modulated bits OVER DISTANCE with noise, interference and attenuation");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. THROUGHPUT
framesFallidos = 0
framesTotales = length(mensaje)
for i = 1:length(mensaje)
    if strcmp(mensaje(i), mensaje1(i)) == 0
        framesFallidos = framesFallidos + 1
    end
end
framesExitosos = framesTotales - framesFallidos
tiempo = (length(t2)-1) / 11000
throuhput = framesExitosos/tiempo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3. ERROR POSIBILITY
FER = framesFallidos/framesTotales

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4. OUTTAGE POSIBILITY
OUTP = (framesExitosos-framesTotales)/framesTotales

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%5. AGE OF INFORMATION
AOI = ((tiempo/framesTotales)*2 - tiempo/framesTotales)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CRC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aux = str2num(str(1,:))

for i = 1:length(str)
    listaSE(i) = bin2dec(str(i,:))
end
    listaCE = mensaje 

errores = 0
for i = 1:length(str) 
    message = listaSE(i)
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
    remainder = bitshift(listaCE(i),divisorDegree);
    remainder = bitor(remainder,CRC_check_value);
    % remainder = bitset(remainder,[2 3]);
    % remainder = bitset(remainder,6);
    
    for k = 1:messageLength
        if bitget(remainder,messageLength+divisorDegree)
            remainder = bitxor(remainder,divisor);
        end
        remainder = bitshift(remainder,1);
    end
    if remainder == 0
        disp('Message is error free.')
    
    else
        disp('Message contains errors.')
        errores = errores + 1
    end
end 