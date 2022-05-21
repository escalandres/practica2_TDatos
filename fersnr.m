clear all
close all
clc

%%%%%Grafica FER SNR%%%%
FER2 = [11.77 3.21 0.76 0.15 0.15];
SNR2 = [1 5 10 20 100];

figure(1)
plot(SNR2,FER2,'LineWidth',2)
xlabel('Signal Noise Rate SNR [dBs]')
ylabel('Frame Error Rate FER [%]')
hold on
grid on
title('FER vs SNR')

