HDR=ScouseTom_getHDR('./nerve030 - Good/N030_EIT001.bdf');

V=sread(HDR,inf,0);
Fs=HDR.SampleRate;
%%
ref_chn= 19; % from pngs plotted
Bad_channels = [1:2,13:29];

% re reference data
Data=V-V(:,ref_chn);

%% James code
Chns_used=1:size(Data,2);
[pks,T_trig]=findpeaks(-V(:,end),'MinPeakHeight',2e5);

% fix the triggers - there are gaps fucking it up

%find where the gap is bigger than the correct value
Sw_in=diff(T_trig);
Sw_True=mode(Sw_in);

%find where the gaps are
Sw_corrected=ceil(Sw_in./Sw_True);
gaps= Sw_corrected > 1;
gapind= find(gaps ==1);


%plug the gaps! this is cheating but who cares/will ever look at this
for gg=1:length(gapind)
    
    pulses_missing=Sw_corrected(gapind(gg)) -1;
    
    new_pulse_vec=(1:pulses_missing)*Sw_True+T_trig(gapind(gg));
    
    T_trig=[T_trig; new_pulse_vec'];
end
T_trig=sort(T_trig);


tau_max=250; % specify in ms

% find max timing between stims
Tmax=ceil(mean(ceil((diff(T_trig)*1000)))/Fs);

tau = min([ tau_max Tmax]);

size_bin=floor(tau*Fs/1000); % convert to the number of samples this is equivalent to


%% spliting into each stim event
cd
% make a Time vector - for plotting
T = [1:size_bin]*1000/Fs;
T = T - T(length(T)/2);

% pre allocate the data

Data_seg=zeros(length(T_trig),size_bin,length(Chns_used));

% loop through every stime event and take the data either side of stim
for iTrig=1:length(T_trig)
    Data_seg(iTrig,:,:)= Data([T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1],:);
end
%%
% for iTrig = 1:length(T_trig)
%
% plot(squeeze(Data_seg(iTrig,:,3)'))
% drawnow
% pause
%
% end

est_seg=Data(floor(T_trig(1)*1.01):ceil(T_trig(2)*0.99),:);

[~,idx]=sort(rms(est_seg),2,'descend');

inj_chn=idx(1:2);
Fc_est=ScouseTom_data_GetCarrier(est_seg(:,inj_chn(1)),Fs);
Fc=round(Fc_est);
inj_chn=sort(inj_chn);

fprintf('Found EIT injection channels %s and %s\n',(HDR.Label{inj_chn(1)}),(HDR.Label{inj_chn(2)}));

%% Kirills bit

good_chn =4; % channel we want to plot 4 is the best one 

%good_trigs=1:140; % which trigger 


Y=Data_seg(1:end,:,good_chn)'; 
xlims = [-10 30];

if mod(size(Y,2),2)==1
    Y = Y(:,2:end);
end

A = detrend(Y(:,1:2:end-2),'constant');
B = detrend(Y(:,2:2:end-2),'constant');
C = detrend(Y(:,3:2:end),'constant');

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
EP = filtfilt(b,a,EP);
EPm= mean(EP,2);

figure
hold on
plot(T,EP,'color',[0.7 0.7 0.7])
plot(T,EPm)
title('EP');
hold off
xlim(xlims)
ylim([-2000,10000])

%%
dV_sig_orig =(A-2*B+C)/4; % kirills linear fit way

BW =125;

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

dV_sigF=filtfilt(FiltBP,dV_sig);

dVdemod=abs(hilbert(dV_sigF));

BV= mean(dVdemod(T > - 80 & T < -20,:));

dV=dVdemod-BV;
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

figure
hold on
plot(T,(dVp),'color',[0.7 0.7 0.7])
plot(T,dVpm)
title('dZp');
hold off
xlim(xlims)
% ylim([-2000,10000])
ylim([-0.1 0.1])





