function check_EP_artefact(EP,Fc)
%

Fs = 16384;
T = 1e3*((1:Fs/2)/Fs-.15);
C = 1e4*sin(Fc*2*pi*T/1e3)';
tClip = 700;
thres = 100;

if mod(size(EP,2),2)==1
    EP = EP(:,2:end);
end

Y = detrend(EP,'constant')+repmat([C,-C],1,size(EP,2)/2);

A = detrend(Y(:,1:2:end-2),'constant');
B = detrend(Y(:,2:2:end-2),'constant');
C = detrend(Y(:,3:2:end),'constant');

% outlier rejection
rY = range(Y);
ff_A = find(abs(range(A)-median(rY))>3*std(rY));
ff_B = find(abs(range(B)-median(rY))>3*std(rY));
ff_C = find(abs(range(C)-median(rY))>3*std(rY));
ff = unique([ff_A,ff_B,ff_C]);
if ~isempty(ff)
    disp(sprintf('Removed %i outliers out of %i',...
                 length(ff),size(A,2)))
    A(:,ff) = [];
    B(:,ff) = [];
    C(:,ff) = [];
end

EP = detrend(-(A+B)/2,'constant');
EP(T>-1.5 & T<1.5,:) = 0;

dV = hilbert_demod(A-2*B+C,Fc,Fs);
mu = mean(mean(dV));
dV = detrend(dV,'constant');
dV(tClip+1:end-tClip,:) = detrend(dV(tClip+1:end-tClip,:));
dV(T>-1.5 & T<1.5,:) = 0;
s = all(abs(dV(T>-10 & T<50,:))<thres);

disp(sprintf('Removed %i noisy trials (threshold = %i uV)',nnz(~s),thres))

figure('Position',[100,0,560,950],'PaperPositionMode','auto')
subplot(2,1,1)
plot(T,EP,'color',.6*[1,1,1]), hold on
plot(T,mean(EP,2),'k','LineWidth',2)
ylim([-1e3,1e3])
xlim([-10,30])
grid on
xlabel('Time (ms)')
ylabel('EP (uV)')

subplot(2,1,2)
plot(T,100*dV(:,s)/mu,'color',.6*[1,1,1]), hold on
plot(T,100*mean(dV(:,s),2)/mu,'k','LineWidth',2)
xlim([-10,30])
ylim([-.1,.1])
grid on
xlabel('Time (ms)')
ylabel('dZ (%)')


