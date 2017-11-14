% clear all

%% CHANGE HERE MARTIN!!!

% change file names
%remove bad repeats, or make it the right number - the variable Y needs to
%be 58
bad_trigs=1:4; % if it breaks make this line 1:3

% 225 hz
% HDR=ScouseTom_getHDR('E:\crabtestingsumsub\170821_07_C31_No1_N1-L1_07_Hooks_1500uA_250uS_EIT_Sum-Sub_Test6.eeg');

% 225 hz longer
% HDR=ScouseTom_getHDR('E:\crabtestingsumsub\170821_13_C31_No2_N2-R1_06_Hooks_1500uA_250uS_EIT_Sum-Sub_Test3.eeg');

% 6k control
% HDR=ScouseTom_getHDR('E:\crabtestingsumsub\170821_22_C31_No2_N2-R1_15_Hooks_1500uA_250uS_EIT_Sum-Sub_Test12.eeg');

% 6k working
% HDR=ScouseTom_getHDR('E:\crabtestingsumsub\170821_21_C31_No2_N2-R1_14_Hooks_1500uA_250uS_EIT_Sum-Sub_Test11.eeg');

% 225 longer again ref 21
% HDR=ScouseTom_getHDR('E:\crabtestingsumsub\170821_11_C31_No2_N2-R1_04_Hooks_1500uA_250uS_EIT_Sum-Sub_Test1.eeg');

% 225 longer again ref 21
% HDR=ScouseTom_getHDR('E:\crabtestingsumsub\170821_12_C31_No2_N2-R1_05_Hooks_1500uA_250uS_EIT_Sum-Sub_Test2.eeg');

% 1k
% HDR=ScouseTom_getHDR('E:\crabtestingsumsub\170821_16_C31_No2_N2-R1_09_Hooks_1500uA_250uS_EIT_Sum-Sub_Test6.eeg');
%
%%

HDR=ScouseTom_getHDR;

Fs = HDR.SampleRate;
[~, fname]=fileparts(HDR.FileName);

%% channel map

switch HDR.TYPE
    case 'BDF' % biosemi file
        map_ = [ 7 8 9:18 25];
        Trigger=ScouseTom_TrigReadChn(HDR);
        TT=ScouseTom_TrigProcess(Trigger,HDR);
    case 'BrainVision'
        
        Trigger = ScouseTom_TrigReadChn(HDR);
        TT=ScouseTom_TrigProcess(Trigger,HDR);
        
    otherwise
        error('Bad HDR');
end

%% read the data
% get the voltages on the desired channels out of the EEG structure
Data = sread(HDR,inf,0);
Chns_used=1:size(Data,2);

bath_chn=1; % free elec channel (if any) used in plotting but nothing else

% Data = Data(:,Good_ch);

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
T = [1:size_bin]*1000/Fs;
T = T - T(length(T)/2);

% pre allocate the data

Data_seg=zeros(length(T_trig),size_bin,length(Chns_used));

% loop through every stime event and take the data either side of stim
for iTrig=1:length(T_trig)
    Data_seg(iTrig,:,:)= Data([T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1],:);
end

xlims = [-10 40];

%% plot all segments

% for iTrig = 1:length(T_trig)
%
% plot(squeeze(Data_seg(iTrig,:,3)'))
% drawnow
% pause
%
% end

%% filtering for EIT Carrier Frequency...

% Fc is the carrier frequency of EIT, BW = ?????

%This line of code detects teh estimated EIT carrier frequency (from the %oscillations within the trace...)

est_seg = Data(floor(T_trig(1)*1.01):ceil(T_trig(2)*0.99),:);

[~,idx] = sort(rms(est_seg),2,'descend');

idx (idx == bath_chn) =[]; % remove bath channel from RMS estimate

inj_chn = idx(1:2);
Fc_est=ScouseTom_data_GetCarrier(est_seg(:,inj_chn(1)),Fs);

inj_chn=sort(inj_chn);

fprintf('Found EIT injection channels %d and %d\n',str2double(HDR.Label{inj_chn(1)}),str2double(HDR.Label{inj_chn(2)}));

good_chn=1:min(inj_chn)-1;


good_chn(good_chn == bath_chn) =[]; % channels which should not artefact - before inj elec

bad_chn = setdiff(Chns_used,[bath_chn good_chn inj_chn]);

% Fc=6000;

Fc = round(Fc_est);

BW =125;
%
% switch Fc
%     case 225
%         BW =125;
%     case 6000
%         BW = 1000;
%     case 2000
%         BW = 500;
% end

DataF=nan(size(Data));
Data_demod=nan(size(Data));


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


N=500;
FiltEP = designfilt('lowpassfir', ...       % Response type
    'FilterOrder',N, ...            % Filter order
    'CutoffFrequency',200, ...    % Frequency constraints
    'DesignMethod','window', ...         % Design method
    'Window','blackmanharris', ...         % Design method options
    'SampleRate',Fs);               % Sample rate

%%


for iChn=1:length(Chns_used)
    
    cur_chn=Chns_used(iChn);
    
    cur_chn_label=str2double(HDR.Label{cur_chn});
    fprintf('Processing #%d elec %d\n',cur_chn,cur_chn_label)
    
    Y=squeeze(Data_seg(:,:,cur_chn))';
    
    Y(:,bad_trigs)=[];
    
    A = detrend(Y(:,1:2:end-2),'constant');
    B = detrend(Y(:,2:2:end-2),'constant');
    C = detrend(Y(:,3:2:end),'constant');
    %% get EP signal
    
    stim_window_ep=[24860 25050];
    
    % summation subtraction to get EPs
    EP = detrend(-(A+B)/2,'constant');
    
    EP=filtfilt(FiltEP,EP);
    
    
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
        try
            cur_dV_sig=(dV_sig(:,iRep));
            [pks,locs]=findpeaks(cur_dV_sig,'MinPeakProminence',max(cur_dV_sig(1:floor(stim_window_dz(1)/2))*0.98));
            
            loc_stim_idx=find(locs > max(stim_window_dz),1);
            loc_good_idx=find(locs < min(stim_window_dz),1,'last');
            
            patch_window=stim_window_dz-(locs(loc_stim_idx)-locs(loc_good_idx));
            
            % figure;
            % hold on
            % plot(cur_dV_sig(stim_window_dz(1):stim_window_dz(2)));
            % plot(cur_dV_sig(patch_window(1):patch_window(2)));
            % hold off
            
            cur_dV_sig(stim_window_dz(1):stim_window_dz(2)) = cur_dV_sig(patch_window(1):patch_window(2));
            dV_sig(:,iRep)=cur_dV_sig;
        catch
            disp('elec fucked up');
        end
    end
    
    %% demodulate to get dZ signal
    
    dV_sigF=dV_sig;
    
    FiltReps=1; % apply filter a different number of times - reduce filter ripple?
    for iFilt = 1:FiltReps
        
        dV_sigF=filtfilt(FiltBP,dV_sigF);
        
    end
    
    dV=abs(hilbert(dV_sigF));
    
    % figure;
    % hold on
    % plot(dV_sigF,'r')
    % plot(dV,'b');
    % hold off
    
    % replace ends of signal to make the next detrends work better
    %
    % trim=700;
    % dV([1:trim],:)=dV([trim+1:2*trim],:);
    % dV([end-trim:end],:)=dV([end-2*trim:end-trim],:);
    
    dV=dV-mean(dV(T > - 100 & T < -2,:));
    
    dV(stim_window_dz(1):stim_window_dz(2),:)=0;
    
    dVm=mean(dV,2);
    %%
    figure;
    
    subplot(2,1,1);
    hold on
    plot(T,dV)
    plot(T,dVm,'k','linewidth',5)
    title(['dZ in channel label: ' num2str(cur_chn_label)]);
    xlabel('Time ms');
    ylabel('Voltage uv');
    xlim(xlims)
    
    y_max = round(max(max((100+dV( T > 5 & T <50,:)))),-2);
    y_min = -round(max(max(100-(-dV( T > 5 & T <50,:)))),-2);
    ylim([y_min y_max])
    
    % tickdiff=mode(diff(get(gca,'XTick')));
    %     extents=max([max(Vsim_full(keep_idx)) Vexp_full(keep_idx)']);
    %     extents=ceil((extents/tickdiff))*tickdiff;
    
    
    subplot(2,1,2);
    hold on
    plot(T,EP)
    plot(T,EPm,'k','linewidth',5)
    title(['EP in channel label: ' num2str(cur_chn_label)]);
    xlabel('Time ms');
    ylabel('Voltage uv');
    xlim(xlims)
    
    
    ep_y_max = round(max(max((100+EP( T > 5 & T <50,:)))),-2);
    ep_y_min = -round(max(max(100-(-EP( T > 5 & T <50,:)))),-2);
    ylim([ep_y_min ep_y_max])
    
    
    drawnow
    
    %%
    
    Out(iChn).EP=EP;
    Out(iChn).EPm=EPm;
    Out(iChn).dV=dV;
    Out(iChn).dVm=dVm;
    Out(iChn).T=T;
    
end


%% Plot all channels


dV_all=nan(size_bin,length(Chns_used));
EP_all=dV_all;

for iChn = 1:length(Chns_used)
    
    dV_all(:,iChn)=Out(iChn).dVm;
    EP_all(:,iChn)=Out(iChn).EPm;
    
end

figure;

subplot(2,1,1);
hold on

title(sprintf('Inj elecs %d %d, @ %d hz and %d uA',str2double(HDR.Label{inj_chn(1)}),str2double(HDR.Label{inj_chn(2)}),Fc,50))
plot(T,EP_all(:,bath_chn),'linewidth',2)
plot(T,EP_all(:,good_chn),'linewidth',3)
plot(T,EP_all(:,inj_chn),'linewidth',0.5)
plot(T,EP_all(:,bad_chn),'linewidth',1)

xlabel('Time ms');
ylabel('EP uv');
xlim(xlims)

ep_y_max = round(max(max((100+EP_all( T > 3 & T <50,:)))),-2);
ep_y_min = -round(max(max((100-EP_all( T > 3 & T <50,:)))),-2);
ylim([ep_y_min ep_y_max])
% ylim([-4100 9000])

legend(HDR.Label)

subplot(2,1,2);
hold on
plot(T,dV_all(:,bath_chn),'linewidth',2)
plot(T,dV_all(:,good_chn),'linewidth',3)
plot(T,dV_all(:,inj_chn),'linewidth',0.5)
plot(T,dV_all(:,bad_chn),'linewidth',1)

xlabel('Time ms');
ylabel('dZ uv');
xlim(xlims)

y_max = round(max(max((10+dV_all( T > 5 & T <50,:)))),-1);
y_min = -round(max(max((10-dV_all( T > 5 & T <50,:)))),-1);
ylim([y_min y_max])
% ylim([-190 440])

legend(HDR.Label)

title(fname,'interpreter','none')













