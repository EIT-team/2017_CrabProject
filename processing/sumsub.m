Fs = 100000;

Y = EPstruc(2).EP(4:end,:,8)';

if mod(size(Y,2),2)==1
    Y = Y(:,2:end);
end

A = detrend(Y(:,1:2:end-2),'constant');
B = detrend(Y(:,2:2:end-2),'constant');
C = detrend(Y(:,3:2:end),'constant');


%% EP from sum sub

%[b,a] = butter(5,200*2/Fs,'low');
%EP = detrend(-(A+B)/2,'constant');
% EP(T>-1.5 & T<1.5,:) = 0;
%EP = filtfilt(b,a,EP);
%EPm= mean(EP,2);


%% EIT from sum sub
dV_sig_orig =(A-2*B+C)/4; % kirills line ar fit way


dV_sigF=dV_sig_orig;

BW = 100;
Fc = 225;
N=1000;
F6dB1=Fc-BW;
F6dB2=Fc+BW;
FiltBP = designfilt('bandpassfir', ...       % Response type
    'FilterOrder',N, ...            % Filter order
    'CutoffFrequency1',F6dB1, ...    % Frequency constraints
    'CutoffFrequency2',F6dB2, ...
    'DesignMethod','window', ...         % Design method
    'Window','blackmanharris', ...         % Design method options
    'SampleRate',Fs);               % Sample rate
    
dV_sigF=filtfilt(FiltBP,dV_sigF);

xmean = mean(dV_sigF);
xcentered = bsxfun(@minus,dV_sigF,xmean);
xampl = abs(hilbert(xcentered));
yupper = bsxfun(@plus,xmean,xampl);

dVdemod=abs(hilbert(dV_sigF));

[hupper,hlower] = envelope(dV_sigF);

plot(mean(hupper,2));
