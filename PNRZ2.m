function [salida] = PNRZ2(h)
    %Codigo de linea: Polar NRZ-L
    %Example:

    %h=[1 0 0 1 1 0 1 0 1 0];

    %PNRZ(h)

    salida = zeros([1 4561556]);

clf;

n=1;

l=length(h);

h(l+1)=1;

lim1=1;
lim2=1001;

ycount = 1;

while n<=length(h)-1;

        t=n-1:0.001:n;
        %t=n-1:n;
    if h(n) == 0

        if h(n+1)==0 

            y=-(t<n)-(t==n);

        else

            y=-(t<n)+(t==n);

        end
        for k = lim1:lim2-1
            salida(k) = salida(k)+y(ycount);
            ycount = ycount+1;
        end
        lim1 = lim1+1001;
        lim2 = lim2+1001;
        ycount = 1;
        %salida(contador)=y;
        d=plot(t,y);grid on;

        title('Line code POLAR NRZ');

        set(d,'LineWidth',2.5);

        hold on;

        axis([0 length(h)-1 -1.5 1.5]);

        %disp('zero');
        %salida(contador)=y;

    else

        if h(n+1)==0

            y=(t<n)-1*(t==n);

        else

            y=(t<n)+1*(t==n);

        end
        %salida(contador)=y;
        for k = lim1:lim2-1
            salida(k) = salida(k)+y(ycount);
            ycount = ycount+1;
        end
        lim1 = lim1+1001;
        lim2 = lim2+1001;
        ycount = 1;
        d=plot(t,y);grid on;

        title('Line code POLAR NRZ');

        set(d,'LineWidth',2.5);

        hold on;

        axis([0 length(h)-1 -1.5 1.5]);

        %disp('one');
        %salida(contador)=y;
    end

    n=n+1;
    contador=contador+1;
    %pause;

end