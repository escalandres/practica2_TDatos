
close all
clear all
clc

contador = 0;
TRecuperado = round(859.5519481);
for k = 249:TRecuperado:132371
    contador = contador+1;
    disp(k)
end

cuca = [65 66 67 68 69 70 71 72];
kuka = native2unicode(cuca);
disp(kuka)

text = 'A long, long time ago.';
binary = dec2bin(text);
str = bin2dec(binary);

textR = native2unicode(str);
disp(textR.')

pA = [1 0 0 0 0 0 1; 1 0 0 0 0 1 0; 1 0 0 0 0 1 1; 1 0 0 0 1 0 0];
pA = string(pA);
pru = zeros([4 1]);

for k=1:4
    pru(k) = bin2dec(pA(k,1)+pA(k,2)+pA(k,3)+pA(k,4)+pA(k,5)+pA(k,6)+pA(k,7));
end
textA = native2unicode(pru);
disp(textA.')

bitsRecovered = zeros([1 length(H)]);
Rcount = 1;
for k=12:4570
    bitsRecovered(Rcount) = bitsRecuperados(k);
    Rcount = Rcount+1;
end