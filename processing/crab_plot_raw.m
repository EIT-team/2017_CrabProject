
HDR=ScouseTom_getHDR;
[~, fname]=fileparts(HDR.FileName);

Trigger=ScouseTom_TrigReadChn(HDR);
TT=ScouseTom_TrigProcess(Trigger,HDR);
Fs=HDR.SampleRate;
%% read the data
% get the voltages on the desired channels out of the EEG structure
Data= sread(HDR,inf,0);

bath_chn=1;
good_chn =4; % channel we want to plot 4 is the best one
other_chn=1:size(Data,2);
other_chn(good_chn)=[];

inj_chn=[3 4];
xlims = [-5 10];

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
T = T - T(length(T)/2);

% pre allocate the data

disp('Segmenting');
Data_seg=zeros(length(T_trig),size_bin,size(Data,2));
% loop through every stime event and take the data either side of stim
for iTrig=1:length(T_trig)
    Data_seg(iTrig,:,:)= detrend(Data((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:));
end

%% Do simple coherent averaging 

% average all the EP chunks across repeats of EPs
Data_avg=detrend(squeeze(mean(Data_seg(start_trial:end,:,:),1)));
%%

figure;

subplot(3,1,1);
hold on

title('All Channel Avgs');
plot(T,Data_avg(:,:)/1000)
xlabel('Time ms');
ylabel('EP mv');
xlim(xlims)

ep_y_max = round(max(max((100+Data_avg( T > 3 & T <xlims(2),good_chn)))),-2);
ep_y_min = -round(max(max((100-Data_avg( T > 3 & T <xlims(2),good_chn)))),-2);
% ylim([ep_y_min ep_y_max])
% ylim([-4100 9000])

% legend(h1,HDR.Label{good_chn})

subplot(3,1,2);
hold on
plot(T,squeeze(Data_seg(:,:,good_chn)/1000),'color',[0.7 0.7 0.7])
plot(T,Data_avg(:,good_chn)/1000)
ylabel('mV')
title(sprintf('EP on elec %s\n',HDR.Label{good_chn}))
hold off
xlim(xlims)

subplot(3,1,3);
hold on
plot(T,Data_avg(:,good_chn)/1000)
plot(T,Data_avg(:,good_chn+1)/1000)
ylabel('mV')
title(sprintf('EP on elec %s\n',HDR.Label{good_chn}))
hold off
xlim(xlims)
% ylim([-1000 1000])

