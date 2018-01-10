
% close all;

%%
HDR=ScouseTom_getHDR;
[~, fname]=fileparts(HDR.FileName);

Trigger=ScouseTom_TrigReadChn(HDR);
TT=ScouseTom_TrigProcess(Trigger,HDR);
Fs=HDR.SampleRate;
%% read the data
% get the voltages on the desired channels out of the EEG structure
Data= sread(HDR,inf,0);

bath_chn=1;
good_chn =3; % channel we want to plot 4 is the best one
other_chn=1:size(Data,2);
other_chn(good_chn)=[];

inj_chn=[3 4];
xlims = [-10 40];

start_trial = 4; % trial we want to start with, for some reason the first few are fucked

%% Triggers

T_trig=TT.Stimulations{1};
% window in ms around event to view
tau_max=250; % specify in ms

% find max timing between stims
Tmax=mean(floor((diff(T_trig)*1000)))/Fs;

tau = min([ tau_max Tmax]);

size_bin=floor(tau*Fs/1000); % convert to the number of samples this is equivalent to


%%

[bep,aep]=butter(3,[125]/(Fs/2),'low');
% DataF_EP=Data;
DataF_EP=filtfilt(bep,aep,Data);



Fc=225;
BW =70;
F6dB1=Fc-BW;
F6dB2=Fc+BW;
% FiltBPq = designfilt('bandpassfir', ...       % Response type
%     'FilterOrder',N, ...            % Filter order
%     'CutoffFrequency1',F6dB1, ...    % Frequency constraints
%     'CutoffFrequency2',F6dB2, ...
%     'DesignMethod','window', ...         % Design method
%     'Window','blackmanharris', ...         % Design method options
%     'SampleRate',Fs);

disp('Filtering for EIT');

[b,a]=butter(3,[F6dB1,F6dB2]/(Fs/2));
DataF_EIT = filtfilt(b,a,Data);
DataF_EIT = abs(hilbert(DataF_EIT));

% [b,a]=butter(1,[45,55]/(Fs/2),'stop');
% DataF_EIT = filtfilt(b,a,DataF_EIT );
%DataF_EIT=ScouseTom_data_DemodHilbert(Data,FiltBPq);




%% spliting into each stim event
% make a Time vector - for plotting
T = (1:size_bin)*1000/Fs;
T = T - T(length(T)/2);

% pre allocate the data

disp('Segmenting');
Data_seg=zeros(floor(length(T_trig)/2)*2,size_bin,size(Data,2));
EITall=Data_seg;
EPall=Data_seg;


% loop through every stime event and take the data either side of stim
for iTrig=1:floor(length(T_trig)/2)*2
    Data_seg(iTrig,:,:)= detrend(Data((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:));
    EITall(iTrig,:,:)= detrend(DataF_EIT((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:));
    EPall(iTrig,:,:)= detrend(DataF_EP((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:));
   % plot(squeeze(Data_seg(iTrig,:,good_chn)));
   % title(num2str(iTrig))
   % pause
end

% figure;
% plot(mean(Data_seg(start_trial:end,:,good_chn),1)); hold on
% plot(detrend(squeeze(mean(Data_seg(start_trial:end,:,4),1))))
%% Do simple coherent averaging 

% average all the EP chunks across repeats of EPs
Data_avg=detrend(squeeze(mean(Data_seg(start_trial:end,:,:),1)));
EIT_avg=detrend(squeeze(mean(EITall(start_trial:end,:,:),1)));
EP_avg=detrend(squeeze(mean(EPall(start_trial:end,:,:),1)));


%%

figure;

subplot(3,1,1);
hold on

title('All Channel Avgs');
hold on
plot(T,squeeze(EPall(start_trial:end,:,good_chn)),'color',[0.7 0.7 0.7])
plot(T,EP_avg(:,good_chn))
ylabel('uV')
title(sprintf('EP on elec %s\n',HDR.Label{good_chn}))
hold off
xlim(xlims)

ep_y_max = round(max(max((100+Data_avg( T > 3 & T <xlims(2),good_chn)))),-2);
ep_y_min = -round(max(max((100-Data_avg( T > 3 & T <xlims(2),good_chn)))),-2);
% ylim([ep_y_min ep_y_max])
% ylim([-4100 9000])

% legend(h1,HDR.Label{good_chn})

subplot(3,1,2);
hold on
plot(T,squeeze(EITall(start_trial:end,:,good_chn)),'color',[0.7 0.7 0.7])
plot(T,EIT_avg(:,good_chn))
ylabel('uV')
title(sprintf('EIT on elec %s\n',HDR.Label{good_chn}))
hold off
xlim(xlims)
%  ylim([-50 50])

subplot(3,1,3);
hold on
plot(T,squeeze(Data_seg(start_trial:end,:,good_chn)),'color',[0.7 0.7 0.7])
plot(T,Data_avg(:,good_chn))
ylabel('uV')
title(sprintf('RAW on elec %s\n',HDR.Label{good_chn}))
hold off
xlim(xlims)


