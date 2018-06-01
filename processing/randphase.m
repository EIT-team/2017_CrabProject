%% read the data
% get the voltages on the desired channels out of the EEG structure
HDR=ScouseTom_getHDR;
[~, fname]=fileparts(HDR.FileName);
Trigger=ScouseTom_TrigReadChn(HDR);
TT=ScouseTom_TrigProcess(Trigger,HDR);
Fs=HDR.SampleRate;
Data = sread(HDR,inf,0) ;
Data(:,18:end) = [] ;

%% Settings - CHANGE THESE FOR YOUR SPECIFIC EXPERIMENTAL PROTOCOL
good_chn = [3 4 5 6 7 8 9 10 11 12] ; 
eit_inj_pairs = [6 7];
eit_freq = 225;
injtime = 72; %seconds
injnum = 5; %injections per second




%% Sort the Stimulation Triggers
j = 1;
k = 1;
g = 0;
for i = 1:length(TT.Stimulations{1})
    trigbois(j,k) = TT.Stimulations{1}(i);
    j = j+1;
    if j > injtime*injnum
        k = k + 1;
        j = 1;
    end
end

%% Other Settings - No need to change
xlims = [-5 40];
other_chn = 1:size(Data,2) ;
start_trial = 4 ; % trial we want to start with, for some reason the first few are fucked

%% Figuring out Timing Windows
T_trig = trigbois(1:injnum*injtime,1); % window in ms around event to view
tau_max = 250; % specify in ms
Tmax = mean(floor((diff(T_trig)*1000)))/Fs; % find max timing between stims
tau = min([ tau_max Tmax]); % choose whichever is smallest
size_bin=floor(tau*Fs/1000); % convert to the number of samples this is equivalent to

%% Estimate the Carrier Frequency
%{
% find first chunk of data
est_trig = (injtime*injnum)/2; % dont use the first one as sometimes its messed up
est_seg=detrend(Data(floor(T_trig(est_trig)*0.99):ceil(T_trig(est_trig+1)*1.01),:));

%find the carrier using this chunk
Fc_est = ScouseTom_data_GetCarrier(est_seg(:,eit_inj_pairs(1)),Fs);
Fc = round(Fc_est);
Fc = 225;
%}

%% Filtering EPs from Raw Data
[bep aep] = butter(3,100/(Fs/2),'low'); % 100 Hz lowpass filter
[bepn aepn] = iirnotch(eit_freq/(Fs/2),(eit_freq/(Fs/2))/35); % EIT frequency specific notch filter
[bepnf aepnf] = iirnotch(50/(Fs/2),(50/(Fs/2))/35); % 50 Hz notch filter
DataF_EPlp = filtfilt(bep,aep,Data); % apply lowpass filter
DataF_EPlpn = filtfilt(bepn,aepn,DataF_EPlp); % apply EIT frequency specific notch filter
DataF_EP = filtfilt(bepnf,aepnf,DataF_EPlpn); % apply 50 Hz notch filter

%{
BW = 100;
N = 2;
F6dB1 = Fc-BW;
F6dB2 = Fc+BW;
FiltBPq = designfilt('bandpassiir', ...       % Response type
    'FilterOrder',N, ...            % Filter order
    'HalfPowerFrequency1',F6dB1, ...    % Frequency constraints
    'HalfPowerFrequency2',F6dB2, ...
    'DesignMethod','Butter', ...         % Design method
    'SampleRate',Fs);

disp('Filtering for EIT');

DataF_EIT=ScouseTom_data_DemodHilbert(Data,FiltBPq);
%}
%% Segmenting Data
T = (1:size_bin)*1000/Fs; % make a time vector
T = T - T(round(length(T)/2));
T_step = T(101)-T(100); % calculate one step size of the time vector
disp('Segmenting');
% preallocate where the segmented data is going
Data_seg=zeros(length(T_trig),size_bin,size(Data,2)); 
EPall=zeros(length(T_trig),size_bin,size(Data,2));
%EITall=zeros(length(T_trig),size_bin,size(Data,2));

% loop through every stime event and take the data either side of stim
for iTrig=1:length(T_trig)
    Data_seg(iTrig,:,:)= Data((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:);
    EPall(iTrig,:,:)= DataF_EP((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:);
   % EITall(iTrig,:,:)= DataF_EIT((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:);
end

%% Average the EPs 

% average all the EP chunks across repeats of EPs
EP_avg=detrend(squeeze(mean(EPall(start_trial:end,:,:),1)));
%EIT_avg=detrend(squeeze(mean(EITall(start_trial:end,:,:),1)));

% find 
%EIT_BV=squeeze(mean(EITall(start_trial:end, T > -100 & T < -30,:),2));
%EIT_BVm=mean(EIT_BV);
%EIT_BVstd=std(EIT_BV);
%% Find the BVs

for i = 1:size(Data_seg,3)
    BV_sig = Data_seg(:,:,i)';
    hupper_bv = abs(hilbert(BV_sig));
    BV_all(:,i) = mean(hupper_bv(T > - 80 & T < -20,:));
    BV_mean = mean(BV_all,1);
    BV_std = std(BV_all);
end

%% EIT Set-up

Y = Data_seg(:,:,eit_inj_pairs(1)-1)'; % choose the data at the artifact free recording electrode
% make the size even if it is not
if mod(size(Y,2),2)==1
    Y = Y(:,2:end);
end
% Band pass filtering
BW = 100; % bandwidth
Fc = eit_freq; % center frequency
F6dB1=Fc-BW; % calculated center frequency minus bandwidth
F6dB2=Fc+BW; % calculated center frequency plus bandwidth
[bbp,abp] = butter(3,[F6dB1 F6dB2]/(Fs/2)); % set up the bandpass filter
Y_bp=filtfilt(bbp,abp,Y); % filter the data with the bandpass

% NOTE: Doing this filtering before the summation subtraction eliminates low
% frequency noise that can throw off the pairing

% set up to match the phases in the randomized summation subtraction
coun = 1;
for i = 2:2:size(Y,2)-1
        fft_s = fft(Y_bp(:,i)); % calculate the fft of the signal
        fft_p = fft(Y_bp(:,i+1)); % calculate the fft of the paired
        [mag_s idx_s] = max(abs(fft_s)); % find where the principle frequency is
        [mag_p idx_p] = max(abs(fft_p)); % find where the principle frequency is
        phase_diff = rad2deg(abs(angle(fft_s(idx_s))-angle(fft_p(idx_p)))); % calculate the phase difference as a check
        % make sure that the phase difference is 180
        if phase_diff <= 181 && phase_diff >=179
            dV_sigF_sumsub(:,coun) = (Y_bp(:,i+1)-Y_bp(:,i))/2;
            ep_sumsub(:,coun) = (Y_bp(:,i+1)+Y_bp(:,i))/2;
            coun = coun + 1;
        end
end

hupper_sumsub = abs(hilbert(dV_sigF_sumsub)); % take the hilbert transform (i.e. envelope) of the paired signals
BV_u_sumsub= mean(hupper_sumsub(T > - 80 & T < -20,:)); % find the mean of the boundary voltages for each pair
dV_sumsub=hupper_sumsub-BV_u_sumsub; % normalize the hilbert transform by the boundary voltage 
dVm_sumsub=mean(dV_sumsub,2); % find the mean normalized hilbert transform
dVp_sumsub=100*(dV_sumsub./BV_u_sumsub); % express the voltage as percentage change
dVpm_sumsub=mean(dVp_sumsub,2); % find the mean percentage change

%{
%% EIT with Random Phase
dV_sig_orig = Y;

figure;
plot(T,Y);
title('Pre Bandpass Signal');
xlim(xlims);

drawnow;

BW = 100;
Fc = eit_freq;
N=1000;
F6dB1=Fc-BW;
F6dB2=Fc+BW;

[bbp,abp] = butter(3,[F6dB1 F6dB2]/(Fs/2));


dV_sigF=filtfilt(bbp,abp,dV_sig_orig);

figure;
plot(T,dV_sigF);
title('Post Bandpass Signal');
xlim(xlims);

drawnow;

hupper = abs(hilbert(dV_sigF));

figure;
plot(T,hupper);
title('Post Hilbert Transform Signal');
xlim(xlims);
drawnow;

BV_u= mean(hupper(T > - 80 & T < -20,:));

dV=hupper-BV_u;
dVm=mean(dV,2);

dVp=100*(dV./BV_u);
dVpm=mean(dVp,2);

figure
hold on
plot(T,(dVp),'color',[0.7 0.7 0.7]);
plot(T,dVpm,'linewidth',3);
title(sprintf('dVp on elec %d\n',eit_inj_pairs(1)-1));
ylabel('%');
xlabel('T ms');
hold off
xlim(xlims);
ylim([-0.4 0.4]);

drawnow;
%{
figure;
colorinc = 0;
hold on
for v = 1:size(dVp,2)
    plot(T,dVp(:,v),'color',[0 0 1-colorinc]);
    colorinc = colorinc + 1/(size(dVp,2));
    
    ylabel('%');
    xlabel('T ms');
    xlim(xlims);
    ylim([-2 2]);
    drawnow;
    pause(0.5);
end
  %}  
%}

%% Plotting Everything
figure;

subplot(2,1,1);
hold on
title(sprintf('Simple Coherent Averaging of EP - Trial: %d\n',i));
h1=plot(T,EP_avg(:,good_chn),'linewidth',3);

xlabel('Time ms');
ylabel('EP uv');
xlim(xlims)

ep_y_max = round(max(max((1000+EP_avg( T > 2 & T <xlims(2),good_chn)))),-2);
ep_y_min = -round(max(max((1000-EP_avg( T > 2 & T <xlims(2),good_chn)))),-2);
ylim([ep_y_min ep_y_max])

legend(h1,HDR.Label{good_chn})

subplot(2,1,2)
errorbar(BV_mean/1000,BV_std/1000);
xlabel('Electrode');
ylabel('BV mV')
title('Boundary Voltages')
drawnow

colorinc = 0;
figure;
hold on
for v = 1:size(EPall,1)
    plot(T,detrend(EPall(v,:,eit_inj_pairs(1)-1)),'color',[1-colorinc 0 0]);
    colorinc = colorinc + 1/(size(EPall,1));
end
xlabel('Time ms');
ylabel('EP uv');
xlim(xlims);
title(sprintf('CAP over all trials on elec %d\n',eit_inj_pairs(1)-1));
hold off
drawnow;

colorinc = 0;
figure;
hold on
for v = 1:size(EPall,1)
    plot(T,detrend(EPall(v,:,eit_inj_pairs(2)+1)),'color',[1-colorinc 0 0]);
    colorinc = colorinc + 1/(size(EPall,1));
end
xlabel('Time ms');
ylabel('EP uv');
xlim(xlims);
title(sprintf('CAP over all trials on elec %d\n',eit_inj_pairs(2)+1));
hold off
drawnow;
    
figure;
plot(T,Y);
title('Pre Bandpass Signal');
xlim(xlims);
drawnow;    
    
figure;
plot(T,dV_sigF_sumsub);
title('Post Bandpass Subtracted Signal');
xlim(xlims);
drawnow;

figure;
plot(T,hupper_sumsub);
title('Post Hilbert Transform Signal');
xlim(xlims);
drawnow;

figure;
hold on
plot(T,(dVp_sumsub),'color',[0.7 0.7 0.7]);
plot(T,dVpm_sumsub,'linewidth',3);
title(sprintf('dVp on elec %d\n',eit_inj_pairs(1)-1));
ylabel('%');
xlabel('T ms');
hold off
xlim(xlims);
drawnow;