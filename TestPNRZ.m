function outCL = PNRZ(h,T)
%Codigo de linea: Polar NRZ-L
%Example:

%h=[1 0 0 1 1 0 1 0 1 0];

%PNRZ(h)

clf;

n=1;

l=length(h);

h(l+1)=1;
outCL = 0;
while n<=length(h)-1;

    t=n-1:T:n;

if h(n) == 0

    if h(n+1)==0 

        y=-(t<n)-(t==n);
        outCL = horzcat(outCL,y);
    else

        y=-(t<n)+(t==n);
        outCL = horzcat(outCL,y);
    end

    d=plot(t,y);grid on;

    title('Line code POLAR NRZ');

    set(d,'LineWidth',2.5);

    hold on;

    axis([0 length(h)-1 -1.5 1.5]);

    %disp('zero');

else

    if h(n+1)==0

        y=(t<n)-1*(t==n);
        outCL = horzcat(outCL,y);
    else

        y=(t<n)+1*(t==n);
        outCL = horzcat(outCL,y);
    end

    d=plot(t,y);grid on;

    title('Line code POLAR NRZ');

    set(d,'LineWidth',2.5);

    hold on;

    axis([0 length(h)-1 -1.5 1.5]);

    %disp('one');

end

n=n+1;

%pause;

end