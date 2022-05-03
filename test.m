close all
clear all
clc

text = 'A long, long time ago. I can still remember how that music used to make me smile. And I knew if I had my chance that I could make those people dance. And maybe they would be happy for a while. But February made me shiver. With every paper I would deliver. Bad news on the doorstep. I could not take one more step. I cannot remember if I cried. When I read about his widowed bride. But something touched me deep inside. The day the music died. So bye-bye, Miss American Pie. Drove my Chevy to the levee, but the levee was dry. And them good old boys were drinking whiskey and rye. Singin This will be the day that I die. This will be the day that I die';
%text = 'A long, long time ago.';
disp('text to binary conversion started')
binary = dec2bin(text);
limit = length(text)*7;
%Vector de 1x4557, donde se guardaran todos los datos binarios de binary
H = bitStreameros ([1 limit]);
count = 1;
for i = 1:length(text)
    for j = 1:7
%         h(count) = str2num(binary(i,j));
        H(count) = str2double(binary(i,j));
        count=count+1;
    end
end
disp('text to binary conversion finished')

Fbit = 100;
Tbit = 1/Fbit;
% k = (length(H)-2)/Fbit;
% n = 0:Tbit:length(H)+k;
disp('Polar NRbitStream codification started')
bitStream=TestPNRbitStream(H,Tbit);
bitStream(1)=[];
bitStream(length(bitStream))=[];
disp('Polar NRbitStream codification finished')

% figure(3)
% plot(n,bitStream,'LineWidth',2.5);grid on;
% axis([-0.1 length(H) -1-0.1 1+0.1])
% title('Line code POLAR NRbitStream');

A = 20;
for k=1:length(bitStream)
    if bitStream(k)==1
        bitStream(k)=A;
    else
        bitStream(k)=-A;
    end
end

% figure(4)
% 
% plot(n,bitStream,'LineWidth',2.5);grid on;
% axis([-0.1 length(H) -A-4 A+4])
% title("Bit stream with bit 1 = "+A+" and bit 0 = "+(-A));

fm = 44100; % frecuencia de muestreo [HertbitStream = muestras por seg];
Tm = 1/fm;
fc = 5010; %HbitStream
tmax=(3*460256)/132301; %tiempo en s
t1 = 0: Tm :tmax;
t1(length(t1)+1)=10.4366;
t1(length(t1)+1)=10.4366;
t1(length(t1)+1)=10.4366;

figure(5)
plot(t1,bitStream)