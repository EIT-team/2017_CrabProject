
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

good_chn = 4; % channel we want to plot 4 is the best one
other_chn=1:size(Data,2);
other_chn(good_chn)=[];

inj_chn=[5 6];
eit_rec = 4;
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

%% spliting into each stim event
% make a Time vector - for plotting
T = (1:size_bin)*1000/Fs;
T = T - T(round(length(T)/2));

%%
disp('Segmenting');
Data_seg=zeros(floor(length(T_trig)/2)*2,size_bin,size(Data,2));
EITall=zeros(floor(length(T_trig)/2)*2,size_bin,size(Data,2));
EPall=zeros(floor(length(T_trig)/2)*2,size_bin,size(Data,2));

%loop through every stime event and take the data either side of stim
for iTrig=1:floor(length(T_trig)/2)*2
    Data_seg(iTrig,:,:)= detrend(Data((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:));
    %EITall(iTrig,:,:)= detrend(DataF_EIT((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:));
    %EPall(iTrig,:,:)= detrend(DataF_EP((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:));
end

[bep,aep]=butter(3,[100]/(Fs/2),'low');
DataF_EP=filtfilt(bep,aep,Data);

%{
Data_clean = zeros(size(Data,1),1);

    p = 2;
    
    trial_max = max(Data(((T_trig(10)-2800):T_trig(10)-800),7));
    trial_min = min(Data(((T_trig(10)-2800):T_trig(10)-800),7));
    
    for g = 1:size(Data,1)-1
       
        delta_V(p) = Data(g+1,7) - Data(g,7);
        cur_V(p) = Data(g,7);
        
        if delta_V(p) > 2000 || cur_V(p) > trial_max || cur_V(p) < trial_min
            Data_clean(g) = cur_V(p-1); 
        end
        if cur_V(p) <= trial_max && cur_V(p) >= trial_min
            Data_clean(g) = cur_V(p);
            p = p + 1;
        end
    end
%}

BW = 100;
Fc = 225;
F6dB1=Fc-BW;
F6dB2=Fc+BW;

disp('Filtering for EIT');

[b,a] = butter(3,[40]/(Fs/2),'low');
[bbp,abp] = butter(3,[F6dB1 F6dB2]/(Fs/2));

DataF_EIT = filtfilt(bbp,abp,Data_seg(:,:,good_chn)');

DataF_dz = abs(hilbert(DataF_EIT));

EITall = filtfilt(b,a,DataF_dz);

for iTrig=1:floor(length(T_trig)/2)*2
    
    EPall(iTrig,:,:)= detrend(DataF_EP((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:));
end


%% Do simple coherent averaging 

% average all the EP chunks across repeats of EPs
Data_avg=detrend(squeeze(mean(Data_seg(start_trial:end,:,:),1)));
EIT_avg=detrend(squeeze(mean(EITall(:,:),2)));
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

subplot(3,1,2);
hold on
plot(T,EITall,'color',[0.7 0.7 0.7])
plot(T,EIT_avg)
ylabel('uV')
title(sprintf('EIT on elec %s\n',HDR.Label{good_chn}))
hold off
xlim(xlims)

subplot(3,1,3);
hold on
plot(T,squeeze(Data_seg(start_trial:end,:,good_chn)),'color',[0.7 0.7 0.7])
plot(T,Data_avg(:,good_chn))
ylabel('uV')
title(sprintf('RAW on elec %s\n',HDR.Label{good_chn}))
hold off
xlim(xlims)


