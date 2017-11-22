
close all;

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
good_chn =[3 4]; % channel we want to plot 4 is the best one
other_chn=1:size(Data,2);
other_chn(good_chn)=[];

inj_chn=[3 4];
xlims = [-20 40];

start_trial = 4; % trial we want to start with, for some reason the first few are fucked

%% Triggers

T_trig=TT.Stimulations{1};
% window in ms around event to view
tau_max=250; % specify in ms

% find max timing between stims
Tmax=mean(floor((diff(T_trig)*1000)))/Fs;

tau = min([ tau_max Tmax]);

size_bin=floor(tau*Fs/1000); % convert to the number of samples this is equivalent to

%% Estimate the Carrier Frequency

% find first chunk of data
est_trig = 5; % dont use the first one as sometimes its messed up
est_seg=detrend(Data(floor(T_trig(est_trig)*1.01):ceil(T_trig(est_trig+1)*0.99),:));

%find the carrier using this chunk
Fc_est=ScouseTom_data_GetCarrier(est_seg(:,inj_chn(1)),Fs);
Fc=round(Fc_est);
% fprintf('Found EIT injection channels %s and %s\n',strtrim(HDR.Label{inj_chn(1)}),strtrim(HDR.Label{inj_chn(2)}));


%% filtering for EP data

N=500;
FcEP=150;

FiltEP = designfilt('lowpassfir', ...       % Response type
    'FilterOrder',N, ...            % Filter order
    'CutoffFrequency',FcEP, ...    % Frequency constraints
    'DesignMethod','window', ...         % Design method
    'Window','blackmanharris', ...         % Design method options
    'SampleRate',Fs);               % Sample rate

disp('Filtering for EPs');

% DataF_EP=filtfilt(FiltEP,Data);
DataF_EP=Data;

BW =125;

N=1000;
F6dB1=Fc-BW;
F6dB2=Fc+BW;
FiltBPq = designfilt('bandpassfir', ...       % Response type
    'FilterOrder',N, ...            % Filter order
    'CutoffFrequency1',F6dB1, ...    % Frequency constraints
    'CutoffFrequency2',F6dB2, ...
    'DesignMethod','window', ...         % Design method
    'Window','blackmanharris', ...         % Design method options
    'SampleRate',Fs);

disp('Filtering for EIT');

DataF_EIT=Data;
% DataF_EIT=ScouseTom_data_DemodHilbert(Data,FiltBPq);



%% spliting into each stim event
% make a Time vector - for plotting
T = (1:size_bin)*1000/Fs;
T = T - T(length(T)/2);

% pre allocate the data

disp('Segmenting');
Data_seg=zeros(length(T_trig),size_bin,size(Data,2));
EPall=Data_seg;
EITall=Data_seg;
% loop through every stime event and take the data either side of stim
for iTrig=1:length(T_trig)
    Data_seg(iTrig,:,:)= Data((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:);
    EPall(iTrig,:,:)= DataF_EP((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:);
    EITall(iTrig,:,:)= DataF_EIT((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:);
end

%% Do simple coherent averaging 

% average all the EP chunks across repeats of EPs
EP_avg=detrend(squeeze(mean(EPall(start_trial:end,:,:),1)));
EIT_avg=detrend(squeeze(mean(EITall(start_trial:end,:,:),1)));

% find 
EIT_BV=squeeze(mean(EITall(start_trial:end, T > -100 & T < -30,:),2));
EIT_BVm=mean(EIT_BV);
EIT_BVstd=std(EIT_BV);

%%

figure;

subplot(3,1,1);
hold on

title('Simple Coherent Averaging');
plot(T,EP_avg(:,other_chn),'color',[0.7 0.7 0.7],'linewidth',1)
h1=plot(T,EP_avg(:,good_chn),'linewidth',3);

xlabel('Time ms');
ylabel('EP uv');
xlim(xlims)

ep_y_max = round(max(max((100+EP_avg( T > 3 & T <xlims(2),good_chn)))),-2);
ep_y_min = -round(max(max((100-EP_avg( T > 3 & T <xlims(2),good_chn)))),-2);
ylim([ep_y_min ep_y_max])
% ylim([-4100 9000])

legend(h1,HDR.Label{good_chn})

subplot(3,1,2);
hold on
plot(T,EIT_avg(:,other_chn),'color',[0.7 0.7 0.7],'linewidth',1)
h2=plot(T,EIT_avg(:,good_chn),'linewidth',3);

xlabel('Time ms');
ylabel('dZ uv');
xlim(xlims)

y_max = round(max(max((10+EIT_avg( T > 5 & T <xlims(2),good_chn)))),-1);
y_min = -round(max(max((10-EIT_avg( T > 5 & T <xlims(2),good_chn)))),-1);
ylim([y_min y_max])
% ylim([-190 440])

legend(HDR.Label)

legend(h2,HDR.Label{good_chn})


subplot(3,1,3)
errorbar(EIT_BVm/1000,EIT_BVstd/1000);
% set(gca,'XTickLabels',strtrim(HDR.Label))
xlabel('Electrode');
ylabel('BV mV')
title('Boundary Voltages')




drawnow
%% Do summation subtraction on good channels

for iChn=1:length(good_chn)
    
    cur_chn=good_chn(iChn);




%good_trigs=1:140; % which trigger
cur_chn_label=str2double(HDR.Label{cur_chn});
fprintf('Summation subtraction on elec %d\n',cur_chn_label)

Y=Data_seg(4:end,:,cur_chn)'; % the first ones are fucked for some reason...

if mod(size(Y,2),2)==1
    Y = Y(:,2:end);
end

A = detrend(Y(:,1:2:end-2),'constant');
B = detrend(Y(:,2:2:end-2),'constant');
C = detrend(Y(:,3:2:end),'constant');
%% plot sanity checks
c=lines(2);


%% EP from sum sub

[b,a] = butter(5,200*2/Fs,'low');
EP = detrend(-(A+B)/2,'constant');
% EP(T>-1.5 & T<1.5,:) = 0;
EP = filtfilt(b,a,EP);
EPm= mean(EP,2);


%% EIT from sum sub
dV_sig_orig =(A-2*B+C)/4; % kirills linear fit way


dV_sigF=dV_sig_orig;

BW =100;

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

FiltReps=1;

for iFilt = 1:FiltReps
    
    dV_sigF=filtfilt(FiltBP,dV_sigF);
    
end

dVdemod=abs(hilbert(dV_sigF));

BV= mean(dVdemod(T > - 80 & T < -20,:));

dV=dVdemod-BV;
% dV(T>-1.5 & T<1.5,:) = 0;
dVm=mean(dV,2);
dVp=100*(dV./BV);
dVpm=mean(dVp,2);


%% plot sum sub

figure

subplot(2,2,1)
hold on
plot(T,A/1000,'color',c(1,:))
plot(T,B/1000,'color',c(2,:))
ylabel('mV')
xlabel('T ms');
xlim(xlims)
title(sprintf('Phase AntiPhase Trials on elec %d\n',cur_chn_label))



subplot(2,2,2)
hold on
plot(T,EP,'color',[0.7 0.7 0.7])
plot(T,EPm)
ylabel('uV')
title(sprintf('EP on elec %d\n',cur_chn_label))
hold off
xlim(xlims)
ylim([-5000,20000])

subplot(2,2,3)
hold on
plot(T,(dV),'color',[0.7 0.7 0.7])
plot(T,dVm)
title(sprintf('dV on elec %d\n',cur_chn_label))
ylabel('uV')
hold off
xlim(xlims)
% ylim([-2000,10000])

subplot(2,2,4)
hold on
plot(T,(dVp),'color',[0.7 0.7 0.7])
plot(T,dVpm)
title(sprintf('dVp on elec %d\n',cur_chn_label))
ylabel('%')
xlabel('T ms')
hold off
xlim(xlims)
drawnow

%% Doing sum sub with patching out stim artefact

dV_sig_patch=dV_sig_orig;



    stim_window_dz=find(T > -1 & T< 1);
    
    % patch bit of the signal which has the stim artefact in it with a "good"
    % sine wave. This reduces the sti artefact in the dZ, but we ignore this
    % anyways
    for iRep = 1:size(dV_sig_patch,2)
        try
            cur_dV_sig=(dV_sig_patch(:,iRep));
            [pks,locs]=findpeaks(cur_dV_sig,'MinPeakProminence',max(cur_dV_sig(1:floor(stim_window_dz(1)/2))*0.98));
            
            loc_stim_idx=find(locs > max(stim_window_dz),1);
            loc_good_idx=find(locs < min(stim_window_dz),1,'last');
            
            patch_window=stim_window_dz-(locs(loc_stim_idx)-locs(loc_good_idx));
            
%             figure;
%             hold on
%             plot(cur_dV_sig(stim_window_dz));
%             plot(cur_dV_sig(patch_window));
%             hold off
            
            cur_dV_sig(stim_window_dz) = cur_dV_sig(patch_window);
            dV_sig_patch(:,iRep)=cur_dV_sig;
        catch
            disp('elec fucked up');
        end
    end

BW =100;

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

FiltReps=1;

dV_sig_patchF=dV_sig_patch;

for iFilt = 1:FiltReps
    
    dV_sig_patchF=filtfilt(FiltBP,dV_sig_patchF);
    
end

dV_patch_demod=abs(hilbert(dV_sig_patchF));

BV_patch= mean(dV_patch_demod(T > - 80 & T < -20,:));

dV_patch=dV_patch_demod-BV_patch;
% dV_patch(T>-1.5 & T<1.5,:) = 0;
dV_patchm=mean(dV_patch,2);
dV_patchp=100*(dV_patch./BV_patch);
dV_patchpm=mean(dV_patchp,2);

%%
figure

subplot(3,1,1)
hold on
plot(T,dV_sig_patch/1000,'color',c(1,:))
plot(T,dV_sig_patchF/1000,'color',c(2,:))
hold off
ylabel('mV')
xlabel('T ms');
xlim(xlims)
title(sprintf('Patched Y %d\n',cur_chn_label))


subplot(3,1,2)
hold on
plot(T,EP,'color',[0.7 0.7 0.7])
plot(T,EPm)
ylabel('uV')
title(sprintf('EP on elec %d\n',cur_chn_label))
hold off
xlim(xlims)
ylim([-5000,20000])

subplot(3,1,3)
hold on
plot(T,(dV_patch),'color',[0.7 0.7 0.7])
plot(T,dV_patchm)
title(sprintf('dV on elec %d\n',cur_chn_label))
ylabel('uV')
hold off
xlim(xlims)
% ylim([-2000,10000])






%% Doing sum sub with trimming the end off

T_trim =1;

dV_sigF=dV_sig_orig(T > T_trim,:);

T2=T(T>T_trim);









%%

    Out(iChn).EP=EP;
    Out(iChn).EPm=EPm;
    Out(iChn).dV=dV;
    Out(iChn).dVm=dVm;
    Out(iChn).dVp=dVp;
    Out(iChn).dVpm=dVpm;
    Out(iChn).T=T;

end

%%

dV_all=nan(size_bin,length(good_chn));
EP_all=dV_all;

for iChn = 1:length(good_chn)
    
    dV_all(:,iChn)=Out(iChn).dVm;
    EP_all(:,iChn)=Out(iChn).EPm;
    
end

figure;
subplot(2,1,1);
plot(T,EP_all)
title('EP')
xlim(xlims)
ylabel('uV')
xlabel('T ms')
legend(HDR.Label{good_chn})
subplot(2,1,2);
plot(T,dV_all)
xlim(xlims)
ylabel('uV')
title('EIT')
xlabel('T ms')
legend(HDR.Label{good_chn})


