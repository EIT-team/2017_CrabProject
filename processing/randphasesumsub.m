%% Read the data
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
injtime = 30; %seconds
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
tau_max = 250; % specify in ms
Tmax = mean(floor((diff(T_trig)*1000)))/Fs; % find max timing between stims
tau = min([ tau_max Tmax]); % choose whichever is smallest
size_bin=floor(tau*Fs/1000); % convert to the number of samples this is equivalent to

%% Filtering EPs from Raw Data
[bep aep] = butter(3,100/(Fs/2),'low'); % 100 Hz lowpass filter
[bepn aepn] = iirnotch(eit_freq/(Fs/2),(eit_freq/(Fs/2))/35); % EIT frequency specific notch filter
[bepnf aepnf] = iirnotch(50/(Fs/2),(50/(Fs/2))/35); % 50 Hz notch filter
DataF_EPlp = filtfilt(bep,aep,Data); % apply lowpass filter
DataF_EPlpn = filtfilt(bepn,aepn,DataF_EPlp); % apply EIT frequency specific notch filter
DataF_EP = filtfilt(bepnf,aepnf,DataF_EPlpn); % apply 50 Hz notch filter

%% Segmenting Data
T = (1:size_bin)*1000/Fs; % make a time vector
T = T - T(round(length(T)/2));
T_step = T(101)-T(100); % calculate one step size of the time vector
disp('Segmenting');
% preallocate where the segmented data is going
Data_seg=zeros(length(T_trig),size_bin,size(Data,2)); 
EPall=zeros(length(T_trig),size_bin,size(Data,2));
% loop through every stime event and take the data either side of stim
for iTrig=1:length(T_trig)
    Data_seg(iTrig,:,:)= Data((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:);
    EPall(iTrig,:,:)= DataF_EP((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:);
end

%% Average the EPs 
% average all the EP chunks across repeats of EPs
EP_avg=detrend(squeeze(mean(EPall(start_trial:end,:,:),1)));

%% Find the BVs
for i = 1:size(Data_seg,3)
    BV_sig = detrend(squeeze(Data_seg(:,:,i)))';
    hupper_bv = abs(hilbert(BV_sig));
    BV_all(:,i) = mean(hupper_bv(T > - 80 & T < -20,:));
    BV_mean = mean(BV_all,1);
    BV_std = std(BV_all);
end

%% EIT Set-up
Y = Data_seg(:,:,eit_inj_pairs(1)-1)'; % choose the data at the artifact free recording electrode
Y_a = Data_seg(:,:,eit_inj_pairs(2)+1)';
% make the size even if it is not
if mod(size(Y,2),2)==1
    Y = Y(:,2:end);
    Y_a = Y_a(:,2:end);
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
            ep_b_eit(:,coun) = (Y(:,i+1)+Y(:,i))/2;
            ep_a_eit(:,coun) = (Y_a(:,i+1)+Y_a(:,i))/2;
            coun = coun + 1;
        end
end

hupper_sumsub = abs(hilbert(dV_sigF_sumsub)); % take the hilbert transform (i.e. envelope) of the paired signals
BV_u_sumsub= mean(hupper_sumsub(T > - 80 & T < -20,:)); % find the mean of the boundary voltages for each pair
dV_sumsub=hupper_sumsub-BV_u_sumsub; % normalize the hilbert transform by the boundary voltage 
dVm_sumsub=mean(dV_sumsub,2); % find the mean normalized hilbert transform
dVp_sumsub=100*(dV_sumsub./BV_u_sumsub); % express the voltage as percentage change
dVpm_sumsub=mean(dVp_sumsub,2); % find the mean percentage change

%% OTHER CALCS
start_ep_b = 5.8544e3;
start_ep_a = 4.0698e3;
start_area_b = 2.6865e4;
start_area_a = 1.9934e4;

ep_b_eit_m = mean(ep_b_eit,2);
ep_a_eit_m = mean(ep_a_eit,2);

ep_b_bl = mean(ep_b_eit_m(T > -20 & T < -10));
ep_a_bl = mean(ep_a_eit_m(T > -20 & T < -10));

ep_b_min = (abs(ep_b_bl - min(ep_b_eit_m(T > 2 & T < 30)))/start_ep_b)*100;
ep_a_min = (abs(ep_a_bl - min(ep_a_eit_m(T > 2 & T < 30)))/start_ep_a)*100;

for i = 1:(15500-12700)
    ep_b_area(i) = T_step*(abs(ep_b_bl - ep_b_eit(12700+i)));
    ep_a_area(i) = T_step*(abs(ep_a_bl - ep_a_eit(12700+i)));
end

ep_b_tot = (sum(ep_b_area)/start_area_b)*100;
ep_a_tot = (sum(ep_a_area)/start_area_a)*100;

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
%{    
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
%}
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