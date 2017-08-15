function x = correct_cos_stimulation(T,x,Fs,Fc);

i_bef=find(T<1.5);
i_edge=max(i_bef);
[val1,point1]=min(abs(x([(i_edge-round(Fs/Fc):i_edge)])));
[val2,point2]=max((x([(i_edge-round(Fs/Fc):i_edge)])));
[val3,point3]=min((x([(i_edge-round(Fs/Fc):i_edge)])));
point1=point1+i_edge-round(Fs/Fc)-1;
point2=point2+i_edge-round(Fs/Fc)-1;
point3=point3+i_edge-round(Fs/Fc)-1;

% mod=find(T<-2);
% 
% L=length (x(mod));
% y=x(mod);
% 
% NFFT = 2^nextpow2(L); % Next power of 2 from length of y
% Y = fft(y,NFFT)/L;
% f = Fs/2*linspace(0,1,NFFT/2+1);
% 
% % Plot single-sided amplitude spectrum.
% [A,k]=max(2*abs(Y(1:NFFT/2+1))); 
% freq=f(k);

A=(val2-val3)/2;
c=(val2+val3)/2;

k1=asin((val1-c)/A);
k2=asin((val2-c)/A);

f=(k1-k2)/(point1-point2);
phi=k2-f*point2;

%y(1:point1)=A*sin(2*pi*1025*[1:point1]/Fs)+c;

mod=find(T<7 & T>2);

L=length (x(mod));
y=x(mod);

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
[Ahh,k]=max(2*abs(Y(1:NFFT/2+1))); 
p=angle(Y(k));
p_samp=round((p+pi)*Fs/(2*pi*Fc));

A=abs(hilbert(x));
A=mean(A(700:(end-700)));
%if (p<0 ) p=p+pi/2; end;
x(1:point1)=A*sin(2*pi*Fc*([1:(point1)])/Fs+p+pi)+c;


% phi=asin(c+val0/A)-2*pi*freq*point1./Fs;
% 
% x(1:point1)=A*sin(2*pi*freq*[1:point1]./Fs+round(phi));
%figure;plot(T,x);