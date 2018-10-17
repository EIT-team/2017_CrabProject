clear all
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
eit_inj_pairs = [4 5];
eit_freq = 225;
eit_bw = 100;
injtime = 20; %seconds
injnum = 1; %injections per second

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
tau_max = 500; % specify in ms
Tmax = mean(floor((diff(T_trig)*1000)))/Fs; % find max timing between stims
tau = min([ tau_max Tmax]); % choose whichever is smallest
size_bin=floor(tau*Fs/1000); % convert to the number of samples this is equivalent to

%% Filtering EPs from Raw Data
%[bep aep] = butter(3,100/(Fs/2),'low'); % 100 Hz lowpass filter
[bepn aepn] = iirnotch(eit_freq/(Fs/2),(eit_freq/(Fs/2))/35); % EIT frequency specific notch filter
[bepnf aepnf] = iirnotch(50/(Fs/2),(50/(Fs/2))/35); % 50 Hz notch filter
%DataF_EPlp = filtfilt(bep,aep,Data); % apply lowpass filter
DataF_EPlpn = filtfilt(bepn,aepn,Data); % apply EIT frequency specific notch filter
DataF_EP = filtfilt(bepnf,aepnf,DataF_EPlpn); % apply 50 Hz notch filter

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

% Band pass filtering
BW = eit_bw; % bandwidth
Fc = eit_freq; % center frequency
F6dB1=Fc-BW; % calculated center frequency minus bandwidth
F6dB2=Fc+BW; % calculated center frequency plus bandwidth
[bbp,abp] = butter(3,[F6dB1 F6dB2]/(Fs/2)); % set up the bandpass filter
Y_bp=filtfilt(bbp,abp,Y); % filter the data with the bandpass

% NOTE: Doing this filtering before the summation subtraction eliminates low
% frequency noise that can throw off the pairing


hupper = abs(hilbert(Y_bp)); % take the hilbert transform (i.e. envelope) of the paired signals
BV_u= mean(hupper(T > - 80 & T < -20,:)); % find the mean of the boundary voltages for each pair
dV=hupper-BV_u; % normalize the hilbert transform by the boundary voltage 
dVm=mean(dV,2); % find the mean normalized hilbert transform
dVp=100*(dV./BV_u); % express the voltage as percentage change
dVpm=mean(dVp,2); % find the mean percentage change

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


