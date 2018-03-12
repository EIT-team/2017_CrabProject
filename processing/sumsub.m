Fs = 100000;
T = (1:25000)*1000/Fs;
T = T - T(length(T)/2);

%Y = EPstruc(2).EP(4:end,:,8)';
Y = EPall(4:end,:,4)';

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

[dVdemod,hlower] = envelope(dV_sigF);

BV= mean(dVdemod(T > - 80 & T < -20,:));

dV=dVdemod-BV;
dVm=mean(dV,2);

figure
hold on
plot(T,(dV),'color',[0.7 0.7 0.7]);
plot(T,dVm,'linewidth',3);
title('dV on elec 4');
ylabel('uV');
xlabel('T ms');
hold off
xlim([-5 30]);
