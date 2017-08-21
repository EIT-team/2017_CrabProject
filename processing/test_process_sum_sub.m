
% clear all
HDR=ScouseTom_getHDR('E:\crabtestingsumsub\170821_07_C31_No1_N1-L1_07_Hooks_1500uA_250uS_EIT_Sum-Sub_Test6.eeg');

%%
Fs = HDR.SampleRate;
[~, fname]=fileparts(HDR.FileName);

Trigger=ScouseTom_TrigReadChn(HDR);
TT=ScouseTom_TrigProcess(Trigger,HDR);

%% channel map

switch HDR.TYPE
    case 'BDF' % biosemi file
        map_ = [ 7 8 9:18 25];
        Trigger=ScouseTom_TrigReadChn(HDR);
        TT=ScouseTom_TrigProcess(Trigger,HDR);
    case 'BrainVision'
        
% %         map_ = [1 2 3 4 5 6 7 8 9 10 11 12];
%         map_ = [1 2 3 4 5 6 7 8 9 10 11 12];
        
        Trigger = ScouseTom_TrigReadChn(HDR);
        TT=ScouseTom_TrigProcess(Trigger,HDR);
        
    otherwise
        error('Bad HDR');
end




%% read the data
% get the voltages on the desired channels out of the EEG structure
Data = sread(HDR,inf,0);
Good_ch=1:size(Data,2);

% Data = Data(:,Good_ch);

%% filtering for EIT Carrier Frequency...

% Fc is the carrier frequency of EIT, BW = ?????

%This line of code detects teh estimated EIT carrier frequency (from the %oscillations within the trace...)
% Fc_est=ScouseTom_data_GetCarrier(Data(:,8),Fs);
Fc=225;

% Fc = round(Fc_est);
BW = 125;

DataF=nan(size(Data));
Data_demod=nan(size(Data));

% bandpass filter -- 1000hz cut off FOR EIT
% [b,a] = butter(5,[Fc-BW Fc+BW ]./(Fs/2),'bandpass');
[b,a]=fir1(500,[Fc-BW Fc+BW ]./(Fs/2),'bandpass');
%DataF = filtfilt(b,a,Data); %comment this to turn it off this only changes eit signal not cap
% Data_demod = abs(hilbert(DataF));

N=1000;
F6dB1=Fc-BW;
F6dB2=Fc+BW;
d = designfilt('bandpassfir', ...       % Response type
    'FilterOrder',N, ...            % Filter order
    'CutoffFrequency1',F6dB1, ...    % Frequency constraints
    'CutoffFrequency2',F6dB2, ...
    'DesignMethod','window', ...         % Design method
    'Window','blackmanharris', ...         % Design method options
    'SampleRate',Fs);               % Sample rate


%low pass filter -- 1000hz cut off FOR REMOVING EIT SIGNAL FROM CAP
[b_ep,a_ep] = butter(5,20000/(Fs/2),'low');
%Data = filtfilt(b_ep,a_ep,Data); %comment this to turn it off



%%
T_trig=TT.Stimulations{1};
% window in ms around event to view
tau_max=500; % specify in ms

% find max timing between stims
Tmax=ceil(mean(ceil((diff(T_trig)*1000)))/Fs);

tau = min([ tau_max Tmax]);

size_bin=floor(tau*Fs/1000); % convert to the number of samples this is equivalent to


%% spliting into each stim event

% make a Time vector - for plotting
T=[1:size_bin]*1000/Fs;
T=T - T(length(T)/2);

% pre allocate the data

Data_seg=zeros(length(T_trig),size_bin,length(Good_ch));

% loop through every stime event and take the data either side of stim
for iTrig=1:length(T_trig)
    Data_seg(iTrig,:,:)= Data([T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1],:);
end


%% plot all segments

% for iTrig = 1:length(T_trig)
%
% plot(squeeze(Data_seg(iTrig,:,3)'))
% drawnow
% pause
%
% end
%% 

cur_chn=1;
cur_chn_label=str2double(HDR.Label{cur_chn});
fprintf('Processing #%d elec %d\n',cur_chn,cur_chn_label)


Y=squeeze(Data_seg(:,:,cur_chn))';

%remove bad repeats, or make it the right number
bad_trigs=1:4;
Y(:,bad_trigs)=[];


A = detrend(Y(:,1:2:end-2),'constant');
B = detrend(Y(:,2:2:end-2),'constant');
C = detrend(Y(:,3:2:end),'constant');
%% get EP signal

stim_window_ep=[24860 25050];

% summation subtraction to get EPs
EP = detrend(-(A+B)/2,'constant');

% replace stim artefact with 0s for ease of plotting
EP(stim_window_ep(1):stim_window_ep(2),:)=0;
EPm=mean(EP,2);

%% get EIT signal 
dV_sig_orig =(A-2*B+C)/4; % kirills linear fit way
% dV_sig =(A-B)/2;

dV_sig=dV_sig_orig;

stim_window_dz=[24860 24980];

% patch bit of the signal which has the stim artefact in it with a "good"
% sine wave. This reduces the sti artefact in the dZ, but we ignore this
% anyways
for iRep = 1:size(dV_sig,2)

cur_dV_sig=(dV_sig(:,iRep));
[pks,locs]=findpeaks(cur_dV_sig,'MinPeakProminence',max(cur_dV_sig(1:floor(stim_window(1)/2))*0.98));

loc_stim_idx=find(locs > max(stim_window),1);
loc_good_idx=find(locs < min(stim_window),1,'last');

patch_window=stim_window_dz-(locs(loc_stim_idx)-locs(loc_good_idx));

% figure;
% hold on
% plot(cur_dV_sig(stim_window_dz(1):stim_window_dz(2)));
% plot(cur_dV_sig(patch_window(1):patch_window(2)));
% hold off

cur_dV_sig(stim_window_dz(1):stim_window_dz(2)) = cur_dV_sig(patch_window(1):patch_window(2));
dV_sig(:,iRep)=cur_dV_sig;

end

%% demodulate to get dZ signal 

dV_sigF=dV_sig;

FiltReps=1; % apply filter a different number of times - reduce filter ripple?
for iFilt = 1:FiltReps
    
    dV_sigF=filtfilt(d,dV_sigF);
    
end

dV=abs(hilbert(dV_sigF));

% figure;
% hold on
% plot(dV_sigF,'r')
% plot(dV,'b');
% hold off

% replace ends of signal to make the next detrends work better

trim=700;
dV([1:trim],:)=dV([trim+1:2*trim],:);
dV([end-trim:end],:)=dV([end-2*trim:end-trim],:);

dV=detrend(dV);

dV(stim_window_dz(1):stim_window_dz(2),:)=0;

dVm=mean(dV,2);
%%
figure;
hold on
plot(T,dV)
plot(T,dVm,'k','linewidth',5)
title(['dZ in channel label: ' num2str(cur_chn_label)]);
xlabel('Time ms');
ylabel('Voltage uv');
xlim([-2 50])

figure;
hold on
plot(T,EP)
plot(T,EPm,'k','linewidth',5)
title(['EP in channel label: ' num2str(cur_chn_label)]);
xlabel('Time ms');
ylabel('Voltage uv');
xlim([-2 50])


% xlim([2.54 2.64]*1e4)































