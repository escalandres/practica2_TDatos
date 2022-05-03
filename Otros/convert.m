function [signal] = convert(h,fc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
lim1=1;
lim2=fc;
tamano=length(h)*fc;
prueba = zeros([1 tamano]);
for i= 1:length(h)
    for j=lim1:lim2
        prueba(j)=prueba(j)+h(i);
    end
    lim1=lim1+fc;
    lim2=lim2+fc;
end
signal=prueba;
end

