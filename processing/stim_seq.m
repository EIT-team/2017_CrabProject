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

Data(:,10:end)=[];

stim_chn = [1 2];
good_chn =[3 4]; % channel we want to plot 4 is the best one
other_chn=1:size(Data,2);
other_chn(good_chn)=[];

inj_chn=[3 4];
xlims = [-20 40];

start_trial = 2; % trial we want to start with, for some reason the first few are fucked

%% Triggers

T_trig=TT.Stimulations{1};
% window in ms around event to view
tau_max=250; % specify in ms

% find max timing between stims
Tmax=mean(floor((diff(T_trig)*1000)))/Fs;

tau = min([ tau_max Tmax]);

size_bin=floor(tau*Fs/1000); % convert to the number of samples this is equivalent to

%% filtering for EP data

N=50;
FcEP=350;

FiltEP = designfilt('lowpassfir', ...       % Response type
    'FilterOrder',N, ...            % Filter order
    'CutoffFrequency',FcEP, ...    % Frequency constraints
    'DesignMethod','window', ...         % Design method
    'Window','blackmanharris', ...         % Design method options
    'SampleRate',Fs);               % Sample rate

disp('Filtering for EPs');

DataF_EP=filtfilt(FiltEP,Data);
% DataF_EP=Data;

%% spliting into each stim event
% make a Time vector - for plotting
T = (1:size_bin)*1000/Fs;
T = T - T(length(T)/2);

% pre allocate the data

disp('Segmenting');
Data_seg=zeros(length(T_trig),size_bin,size(Data,2));
EPall=Data_seg;
% loop through every stime event and take the data either side of stim
for iTrig=1:length(T_trig)
    Data_seg(iTrig,:,:)= Data((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:);
    EPall(iTrig,:,:)= DataF_EP((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:);
end



%% Do simple coherent averaging 

% assume 5 reps and 4 stim levels

% stim level 1

stim_reps=5;
stim_levels=4;


for iStim = 1:stim_levels
start_idx=start_trial+(stim_reps*(iStim-1));
    cur_idx=start_idx:start_idx+stim_reps-1
    
% average all the EP chunks across repeats of EPs
EP_avg(:,:,iStim)=detrend(squeeze(mean(EPall(cur_idx,:,:),1)));

end


%%

figure;

subplot(2,1,1);
hold on

title('Stim Chn');
% plot(T,EP_avg(:,other_chn,:),'color',[0.7 0.7 0.7],'linewidth',1)
h1=plot(T,squeeze(EP_avg(:,stim_chn(1),:)),'linewidth',3);

xlabel('Time ms');
ylabel('Stim uV');
xlim(xlims)

ep_y_max = round(max(max((1000+EP_avg( T > -2 & T <xlims(2),stim_chn(1),:)))),-2);
ep_y_min = -round(max(max((1000-EP_avg( T > -2 & T <xlims(2),stim_chn(1),:)))),-2);
ylim([ep_y_min ep_y_max])
% ylim([-4100 9000])

legend(h1,'Stim1','Stim2','Stim3','Stim4')
hold off


subplot(2,1,2);
hold on

title('EP chn');
% plot(T,EP_avg(:,other_chn,:),'color',[0.7 0.7 0.7],'linewidth',1)
h1=plot(T,squeeze(EP_avg(:,good_chn(1),:)),'linewidth',3);
% h2=plot(T,squeeze(EP_avg(:,good_chn(2),:)),'linewidth',3);

xlabel('Time ms');
ylabel('Stim uV');
xlim(xlims)

ep_y_max = round(max(max((1000+EP_avg( T > -2 & T <xlims(2),good_chn(1),:)))),-2);
ep_y_min = -round(max(max((1000-EP_avg( T > -2 & T <xlims(2),good_chn(1),:)))),-2);
ylim([ep_y_min ep_y_max])
% ylim([-4100 9000])

legend(h1,'Stim1','Stim2','Stim3','Stim4')
hold off



drawnow