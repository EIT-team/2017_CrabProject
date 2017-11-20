HDR=ScouseTom_getHDR('E:\crabcheck\225 Hz\170901_02_C33_No5_N4-R2_02_LinearArray_1500uA_250uS_EIT_50uA_225Hz_1.eeg');

V=sread(HDR,inf,0);
Trigger = ScouseTom_TrigReadChn(HDR);
TT=ScouseTom_TrigProcess(Trigger,HDR);
Fs=HDR.SampleRate;
%%
% ref_chn= 10; % from pngs plotted
% re reference data
% Data=V-V(:,ref_chn);
Data=V;

%% James code
Chns_used=1:size(Data,2);
T_trig=TT.Stimulations{1};

tau_max=250; % specify in ms

% find max timing between stims
Tmax=ceil(mean(ceil((diff(T_trig)*1000)))/Fs);

tau = min([ tau_max Tmax]);

size_bin=floor(tau*Fs/1000); % convert to the number of samples this is equivalent to


%% spliting into each stim event
% make a Time vector - for plotting
T = [1:size_bin]*1000/Fs;
T = T - T(length(T)/2);

% pre allocate the data

Data_seg=zeros(length(T_trig),size_bin,length(Chns_used));

% loop through every stime event and take the data either side of stim
for iTrig=1:length(T_trig)
    Data_seg(iTrig,:,:)= Data((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:);
end
%%
% for iTrig = 1:length(T_trig)
% 
% plot(squeeze(Data_seg(iTrig,:,3)'))
% drawnow
% pause
% 
% end


%%
est_seg=Data(floor(T_trig(1)*1.01):ceil(T_trig(2)*0.99),:);

[~,idx]=sort(rms(est_seg),2,'descend');

inj_chn=idx(1:2);
Fc_est=ScouseTom_data_GetCarrier(est_seg(:,inj_chn(1)),Fs);
Fc=round(Fc_est);
inj_chn=sort(inj_chn);

figure
plot((detrend(est_seg(:,:))))

figure
plot(rms(detrend(est_seg(:,:))))


fprintf('Found EIT injection channels %s and %s\n',strtrim(HDR.Label{inj_chn(1)}),strtrim(HDR.Label{inj_chn(2)}));

%% Kirills bit
xlims = [-30 40];
good_chn =3; % channel we want to plot 4 is the best one

%good_trigs=1:140; % which trigger


Y=Data_seg(4:end,:,good_chn)'; % the first ones are fucked for some reason...


if mod(size(Y,2),2)==1
    Y = Y(:,2:end);
end

figure;plot(T,Y);
xlim(xlims)



A = detrend(Y(:,1:2:end-2),'constant');
B = detrend(Y(:,2:2:end-2),'constant');
C = detrend(Y(:,3:2:end),'constant');
%%
figure
subplot(3,1,1)
plot(T,A)
xlim(xlims)
title('A')
subplot(3,1,2)
plot(T,B)
title('B')
xlim(xlims)
subplot(3,1,3)
plot(T,C)
xlim(xlims)
title('C')

[b,a] = butter(5,200*2/Fs,'low');
EP = detrend(-(A+B)/2,'constant');
EP(T>-1.5 & T<1.5,:) = 0;
% EP = filtfilt(b,a,EP);
EPm= mean(EP,2);

figure
hold on
plot(T,EP,'color',[0.7 0.7 0.7])
plot(T,EPm)
title('EP');
hold off
xlim(xlims)
ylim([-5000,20000])

%%
dV_sig_orig =(A-2*B+C)/4; % kirills linear fit way

dV_sigF=dV_sig_orig;

BW =200;

N=100;
F6dB1=Fc-BW;
F6dB2=Fc+BW;
FiltBP = designfilt('bandpassfir', ...       % Response type
    'FilterOrder',N, ...            % Filter order
    'CutoffFrequency1',F6dB1, ...    % Frequency constraints
    'CutoffFrequency2',F6dB2, ...
    'DesignMethod','window', ...         % Design method
    'Window','blackmanharris', ...         % Design method options
    'SampleRate',Fs);               % Sample rate

FiltReps=10;

for iFilt = 1:FiltReps
    
    dV_sigF=filtfilt(FiltBP,dV_sigF);
    
end

dVdemod=abs(hilbert(dV_sigF));

BV= mean(dVdemod(T > - 80 & T < -20,:));

dV=dVdemod-BV;
dV(T>-1.5 & T<1.5,:) = 0;
dVm=mean(dV,2);
dVp=100*(dV./BV);
dVpm=mean(dVp,2);

figure
hold on
plot(T,(dV),'color',[0.7 0.7 0.7])
plot(T,dVm)
title('dZ');
hold off
xlim(xlims)
% ylim([-2000,10000])

% figure
% hold on
% plot(T,(dVp),'color',[0.7 0.7 0.7])
% plot(T,dVpm)
% title('dZp');
% hold off
% xlim(xlims)
% ylim([-2000,10000])
% ylim([-0.1 0.1])





