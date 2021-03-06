clear all

%% Read the data
% get the voltages on the desired channels out of the EEG structure
% you must have run the load_data_master installer before this part will
% work!
HDR=ScouseTom_getHDR;
[~, fname]=fileparts(HDR.FileName);
Trigger=ScouseTom_TrigReadChn(HDR);
TT=ScouseTom_TrigProcess(Trigger,HDR);
Fs=HDR.SampleRate;
Data = sread(HDR,inf,0);
Data(:,27:end) = []; %you need to change 27 to whatever number of channels you have in the recording

%% Settings - CHANGE THESE FOR YOUR SPECIFIC EXPERIMENTAL PROTOCOL
good_chn = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19] ; %the number of channels that will be used
eit_inj_pairs = [11 12]; %where the current injection channels were in the experiment
injtime = 20; % Total time of recording
injnum = 1; % Number of injections per second
eit_freq = 625;
eit_bw = 500;
eit_cur = 10; %EIT current used - this does not change any of the data processing, only some of the values given for the leakage current
%% Sort the Stimulation Triggers
j = 1;
k = 1;

for i = 1:length(TT.Stimulations{1})
    trigbois(j,k) = TT.Stimulations{1}(i);
    j = j+1;
    if j > injtime*injnum
        k = k + 1;
        j = 1;
    end
end
    
%% Other Settings - No need to change
xlims = [-5 40]; %time in ms - can change this if you want to

%% Figuring out Timing Windows
T_trig = trigbois(1:injnum*injtime,1); % 
tau_max = 500; % window in ms around event to view - specify in ms
Tmax = mean(floor((diff(T_trig)*1000)))/Fs; % find max timing between stims
tau = min([tau_max Tmax]); % choose whichever is smallest
size_bin=floor(tau*Fs/1000); % convert to the number of samples this is equivalent to

%% Filtering EPs from Raw Data
[bepn, aepn] = iirnotch(eit_freq/(Fs/2),1.2857e-04);% EIT frequency specific notch filter
[bepnf aepnf] = iirnotch(50/(Fs/2),(50/(Fs/2))/35); % 50 Hz notch filter

[beplp, aeplp] = butter(3,(3000)/(Fs/2),'low'); % 3000 Hz low pass filter
[bephp, aephp] = butter(3,(10)/(Fs/2),'high'); % 3 Hz High pass filter

DataF_EPlp = filtfilt(beplp, aeplp, detrend(squeeze(Data))); %Lowpass
DataF_EPhp = filtfilt(bephp, aephp, DataF_EPlp); %highpass
DataF_EPeitn = filtfilt(bepn,aepn,DataF_EPhp); % apply EIT frequency specific notch filter
DataF_EP = filtfilt(bepnf,aepnf,DataF_EPeitn); % apply 50 Hz notch filter

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
EP_avg=detrend(squeeze(mean(EPall(2:injtime-1,:,:),1)));

%% Find the BVs
for i = 1:size(Data_seg,3)
    BV_sig = detrend(squeeze(Data_seg(:,:,i)))';
    hupper_bv = abs(hilbert(BV_sig));
    BV_all(:,i) = mean(hupper_bv(T > - 80 & T < -20,:));
    BV_mean = mean(BV_all,1);
    BV_std = std(BV_all);
end

%% EIT Set-up
Y = detrend(squeeze(Data_seg(:,:,eit_inj_pairs(1)-1)')); % choose the data at the artifact free recording electrode
Y_a = detrend(squeeze(Data_seg(:,:,eit_inj_pairs(2)+1)')); %This is the data from the electrode after EIT current, important for some calc

% make the size even if it is not
if mod(size(Y,2),2)==1
    Y = Y(:,2:end);
    Y_a = Y_a(:,2:end);
end

mY = mean(Y((T > - 80 & T < -20),:),1);
Y_bp = Y - mY; % remove DC offset

% Band pass filtering
Fc = eit_freq; % center frequency
BW = eit_bw; % bandwidth
F6dB1=Fc-BW; % calculated center frequency minus bandwidth
F6dB2=Fc+BW; % calculated center frequency plus bandwidth
[bbp,abp] = butter(3,[F6dB1 F6dB2]/(Fs/2)); % set up the bandpass filter
%Y_bp=filtfilt(bbp,abp,Y); % filter the data with the bandpass

mep_a_eit = mean(EPall(:,(T > - 80 & T < -20),eit_inj_pairs(2)+1)'); %calc mean offset of ep after EIT injection
mep_b_eit = mean(EPall(:,(T > - 80 & T < -20),eit_inj_pairs(1)-1)'); %calc mean offset of ep before EIT injection
ep_b_eit = EPall(:,:,eit_inj_pairs(1)-1)' - mep_b_eit; %calc eps before EIT injection
ep_a_eit = EPall(:,:,eit_inj_pairs(2)+1)' - mep_a_eit; % calc eps after EIT injection
mean_ep_b_eit = mean(ep_b_eit,2); % mean ep for before EIT inj
% set up to match the phases in the randomized summation subtraction

if mod(injtime,2) == 0 %if inj time is even needs to start at trial 2 because arduino code is messy
    start_trial = 2;
end
if mod(injtime,2) == 1;%if inj time is odd needs to start at trial 1
    start_trial = 1;
end

coun = 1;% counting variable
for i = start_trial:2:size(Y,2)-1
        fft_s = fft(Y_bp(:,i)); % calculate the fft of the signal
        fft_p = fft(Y_bp(:,i+1)); % calculate the fft of the paired
        [mag_s idx_s] = max(abs(fft_s)); % find where the principle frequency is
        [mag_p idx_p] = max(abs(fft_p)); % find where the principle frequency is
        phase_ang(i) = rad2deg(abs(angle(fft_s(idx_s))));%calc phase angle of first inj
        phase_ang(i+1) = rad2deg(abs(angle(fft_p(idx_p))));%calc phase angle of next inj
        phase_diff(coun) = rad2deg(abs(angle(fft_s(idx_s))-angle(fft_p(idx_p)))); % calculate the phase difference as a check
        
        [r_a, lag_a] = xcorr(ep_a_eit(:,i),ep_a_eit(:,i+1)); %calculation of cross correlation of the ep before and after EIT - not currently used
        [~,I_a] = max(abs(r_a));
        lagDiff_a(coun) = lag_a(I_a)/Fs;
        [r_b, lag_b] = xcorr(ep_b_eit(:,i),ep_b_eit(:,i+1));
        [~,I_b] = max(abs(r_b));
        lagDiff_b(coun) = lag_b(I_b)/Fs;
        
        % make sure that the phase difference is 180
        if phase_diff(coun) <= 190 && phase_diff(coun) >= 160 %check if phase diff isnt near 180
            
            %meowsig = ep_b_eit(:,i)-ep_b_eit(:,i+1);
            %testsig = abs(ep_b_eit(:,i)-ep_b_eit(:,i+1));
            %vdiff_b(:,coun) = abs(ep_b_eit(:,i)-ep_b_eit(:,i+1));
            %threshmax = mean(max(testsig((T > - 80 & T < -20))));
            %testsig(testsig < threshmax) = 0;
            %testsig((T > -124 & T < 2)) = 0;
            %vdiff_b(:,coun) = testsig;
            %vdiff_c(:,coun) = meowsig;
            %dzdiff_b(:,coun) = 100*(vdiff_b(:,coun)/BV_mean(eit_inj_pairs(1)-1));
            
            dV_sigF_sumsub(:,coun) = (Y_bp(:,i+1)-Y_bp(:,i))/2;
            
            coun = coun + 1;
        end
end

%{
hupper_dzdiff = abs(hilbert(dzdiff_b));
BV_u_dzdiff = mean(hupper_dzdiff(T > -80 & T < -20,:));
dzdiff_norm = hupper_dzdiff - BV_u_dzdiff;
BV_pre_dzdiff = mean(dzdiff_b(T > -80 & T < -20,:));
dzdiff_pre_norm = dzdiff_b - BV_pre_dzdiff;
%}

dV_sigF_sumsub_pre=filtfilt(bbp,abp,dV_sigF_sumsub); % filter the data with the bandpass
hupper_sumsub = abs(hilbert(dV_sigF_sumsub_pre)); % take the hilbert transform (i.e. envelope) of the paired signals
BV_u_sumsub= mean(hupper_sumsub(T > - 80 & T < -20,:)); % find the mean of the boundary voltages for each pair
dV_sumsub=hupper_sumsub-BV_u_sumsub; % normalize the hilbert transform by the boundary voltage 
dVm_sumsub=mean(dV_sumsub,2); % find the mean normalized hilbert transform
dVp_sumsub=100*(dV_sumsub./BV_u_sumsub); % express the voltage as percentage change
dVpm_sumsub=mean(dVp_sumsub,2); % find the mean percentage change

EPstruc(1).EP = Data_seg;
EPstruc(1).EP_avg = EP_avg;
EPstruc(1).BVrecmean = BV_mean(eit_inj_pairs(1)-1);
EPstruc(1).BVlasteitmean = BV_mean(eit_inj_pairs(2));
EPstruc(1).EPmaxf = min(EP_avg(T >= 1 & T <= 75,eit_inj_pairs(1)-1));
EPstruc(1).EPmaxeit = min(EP_avg(T >= 1 & T <= 75,eit_inj_pairs(2)));
%EPstruc(1).EPmaxinref = max(EP_avg(T >= 21 & T <= 75,18));
%EPstruc(1).EPcv_four = (4*3) / T(12600+find(EP_avg(12600:20000,4) <= min(EP_avg(12600:20000,4))));
%EPstruc(1).EPcv_sev = (4*6) / T(12600+find(EP_avg(12600:20000,7) <= min(EP_avg(12600:20000,7))));
%EPstruc(1).EPcv_elev = (4*10) / T(12600+find(EP_avg(12600:20000,11) <= min(EP_avg(12600:20000,11))));
t_cent_indx = find(EP_avg(T >= 1 & T <= 75,eit_inj_pairs(2)) <= EPstruc(1).EPmaxeit);
EPstruc(1).EPeiteltime = T(find(T==1)+t_cent_indx);
EPstruc(1).EPmaxftime = T(find(T==1)+find(EP_avg(T >= 1 & T <= 75,eit_inj_pairs(1)-1) <= EPstruc(1).EPmaxf));
EPstruc(1).nervez = (BV_mean(eit_inj_pairs(1))-BV_mean(eit_inj_pairs(2)))/(eit_cur);
EPstruc(1).backcur_fivefour = (BV_mean(eit_inj_pairs(1)-1)-BV_mean(eit_inj_pairs(1)))/(EPstruc(1).nervez);
EPstruc(1).backcur_sixseven = (BV_mean(eit_inj_pairs(2))-BV_mean(eit_inj_pairs(2)+1))/(EPstruc(1).nervez);

pstim_line(size(dVm_sumsub),1) = mean(dVm_sumsub(T >= -5 & T <= -1));
dVa = [];
ep_area = [];
for l=1:4000
    epdelta = mean_ep_b_eit(l+find(T==1))-mean(mean_ep_b_eit(T >= -5 & T <= -1));
    dVdelta = dVm_sumsub(l+find(T==0))-mean(dVm_sumsub(T >= -5 & T <= -1));
    dVa(l) = abs(T_step*dVdelta);
    ep_area(l) = abs(T_step*epdelta);
end

dVstruc(1).dV = dV_sumsub;
dVstruc(1).dVp = dVp_sumsub;
dVstruc(1).dVa = sum(dVa);
dVstruc(1).dVpmmin = min(dVpm_sumsub(find(T==1)+t_cent_indx-200:find(T==1)+t_cent_indx+200)); %finding max dZ within 2 ms window from CAP on last EIT electrode
dVstruc(1).dVpmmax = max(dVpm_sumsub(find(T==1)+t_cent_indx-200:find(T==1)+t_cent_indx+200)); %finding max dZ within 2 ms window from CAP on last EIT electrode
dVstruc(1).dVpmintime = T(find(T==1)+t_cent_indx-200+find(dVpm_sumsub(find(T==1)+t_cent_indx-200:find(T==1)+t_cent_indx+200) <= dVstruc(1).dVpmmin));
dVstruc(1).dVpmaxtime = T(find(T==1)+t_cent_indx-200+find(dVpm_sumsub(find(T==1)+t_cent_indx-200:find(T==1)+t_cent_indx+200) <= dVstruc(1).dVpmmax));
dVstruc(1).dVmmin = min(dVm_sumsub(find(T==1)+t_cent_indx-200:find(T==1)+t_cent_indx+200)); %finding max dZ within 2 ms window from CAP on last EIT electrode
dVstruc(1).dVmmax = max(dVm_sumsub(find(T==1)+t_cent_indx-200:find(T==1)+t_cent_indx+200)); %finding max dZ within 2 ms window from CAP on last EIT electrode
dVstruc(1).dVmintime = T(find(T==1)+t_cent_indx-200+find(dVm_sumsub(find(T==1)+t_cent_indx-200:find(T==1)+t_cent_indx+200) <= dVstruc(1).dVmmin));
dVstruc(1).dVmaxtime = T(find(T==1)+t_cent_indx-200+find(dVm_sumsub(find(T==1)+t_cent_indx-200:find(T==1)+t_cent_indx+200) <= dVstruc(1).dVmmax));

%% OTHER CALCS
%{
start_ep_b = 4.3354e3;
start_ep_a = 2.2883e3;
start_area_b = 6.6108e4;
start_area_a = 1.5299e4;

ep_b_eit_m = mean(ep_b_eit,2);
ep_a_eit_m = mean(ep_a_eit,2);

ep_b_bl = mean(ep_b_eit_m(T > -20 & T < -10));
ep_a_bl = mean(ep_a_eit_m(T > -20 & T < -10));

ep_b_min = (abs(ep_b_bl - min(ep_b_eit_m(T > 2 & T < 30)))/start_ep_b)*100;
ep_a_min = (abs(ep_a_bl - min(ep_a_eit_m(T > 2 & T < 30)))/start_ep_a)*100;

for i = 1:(17500-12700)
    ep_b_area(i) = T_step*(abs(ep_b_bl - ep_b_eit(12700+i)));
    ep_a_area(i) = T_step*(abs(ep_a_bl - ep_a_eit(12700+i)));
end

ep_b_tot = (sum(ep_b_area)/start_area_b)*100;
ep_a_tot = (sum(ep_a_area)/start_area_a)*100;
%}
%% Plotting Everything

figure;
subplot(2,1,1);
set(gca,'FontSize',18)
hold on
title(sprintf('Simple Coherent Averaging of EP - Trial: %d\n',i),'FontSize', 18);
h1=plot(T,(EP_avg(:,good_chn))/1000,'linewidth',4);
xlabel('Time (ms)','FontSize', 18);
ylabel('EP (mV)','FontSize', 18);
xlim(xlims)
ep_y_max = (round(max(max((1000+EP_avg( T > 2 & T <xlims(2),good_chn)))),-2))/1000;
ep_y_min = (-round(max(max((1000-EP_avg( T > 2 & T <xlims(2),good_chn)))),-2))/1000;
ylim([ep_y_min ep_y_max])
legend(h1,HDR.Label{good_chn})
subplot(2,1,2);
set(gca,'FontSize',18)
errorbar(BV_mean/1000,BV_std/1000,'linewidth',3);
xlabel('Electrode','FontSize', 18);
ylabel('BV (mV)','FontSize', 18)
title('Boundary Voltages','FontSize', 18)
drawnow

colorinc = 0;
figure;
set(gca,'FontSize',18)
hold on
for v = 1:size(EPall,1)
    plot(T,(detrend(EPall(v,:,eit_inj_pairs(1)-1)))/1000,'color',[1-colorinc 0 0], 'linewidth',3);
    colorinc = colorinc + 1/(size(EPall,1));
end

%plot(T,detrend(EPall(:,:,3)')/1000,'color',[0.7 0.7 0.7],'linewidth',3);

%plot(T,(EP_avg(:,3))/1000,'linewidth',6);

xlabel('Time (ms)','FontSize', 18);
ylabel('EP (mV)','FontSize', 18);
xlim(xlims);
ylim([-8 8]);
title(sprintf('CAP over all trials on elec %d\n',eit_inj_pairs(1)-1),'FontSize', 18);
hold off
drawnow;

colorinc = 0;
figure;
set(gca,'FontSize',18)
hold on
for v = 1:size(EPall,1)
    plot(T,(detrend(EPall(v,:,eit_inj_pairs(2)+1)))/1000,'color',[1-colorinc 0 0],'linewidth',3);
    colorinc = colorinc + 1/(size(EPall,1));
end
xlabel('Time (ms)','FontSize', 18);
ylabel('EP (mV)','FontSize', 18);
xlim(xlims);
title(sprintf('CAP over all trials on elec %d\n',eit_inj_pairs(2)+1),'FontSize', 18);
hold off
drawnow;

%meowmod = 100*(vdiff_b/BV_mean(eit_inj_pairs(1)-1));
%meowmean = mean(meowmod(T > -80 & T < -20,:));
%meownorm = meowmod - meowmean;

figure;
set(gca,'FontSize',16)
hold on
plot(T,(dVp_sumsub),'color',[0.7 0.7 0.7],'linewidth', 3);
plot(T,dVpm_sumsub,'linewidth',6);

%plot(T,(filtfilt(blp,alp,meownorm)),'color',[0.7 0.7 0.7],'linewidth',3);
%plot(T,mean(filtfilt(blp,alp,meownorm),2),'linewidth',6);
%plot(T,(filtfilt(blp,alp,dzdiff_norm)),'color',[0.7 0.7 0.7],'linewidth',3);
%plot(T,mean(filtfilt(blp,alp,dzdiff_norm),2),'linewidth',6);

title(sprintf('dVp on elec %d\n',eit_inj_pairs(1)-1),'FontSize', 18);
ylabel('%','FontSize', 18);
xlabel('T ms','FontSize', 18);
hold off
xlim(xlims);
drawnow;

%making the variables easy to read off for importing to spreadsheet
%capmod = filtfilt(blp,alp,dzdiff_norm);
%mean_capmod = mean(filtfilt(blp,alp,dzdiff_norm),2);

if eit_freq < 600 || eit_freq > 900 
    intstuff.max_sumsumdzneg = dVstruc(1).dVpmmin;
    intstuff.max_sumsubdzuv = dVstruc(1).dVmmin;
    intstuff.sumsubmintime = dVstruc(1).dVpmintime;
end
if eit_freq > 600 && eit_freq < 900
    intstuff.max_sumsumdzpos = dVstruc(1).dVpmmax;
    intstuff.sumsubmaxtime = dVstruc(1).dVpmaxtime;
end
intstuff.sumsubarea = dVstruc(1).dVa;
%intstuff.max_capmod = max(mean_capmod(T > 3 & T < 30));
%intstuff.capmodtime = T(find(T==1)+t_cent_indx-200+find(mean_capmod(find(T==1)+t_cent_indx-200:find(T==1)+t_cent_indx+200) >= intstuff.max_capmod));
intstuff.sumsubnoise = mean(max(dVp_sumsub(T > -15 & T < -5,:))-min(dVp_sumsub(T > -15 & T < -5,:)));
%intstuff.capmodnoise = std(capmod(T > -20 & T < -5,:));
intstuff.max_capbefore = EPstruc(1).EPmaxf;
intstuff.beicaptime = EPstruc(1).EPmaxftime;
intstuff.max_capafter = EPstruc(1).EPmaxeit;
intstuff.eicaptime = EPstruc(1).EPeiteltime;
intstuff.capnoise = mean(std(ep_b_eit((T > -20 & T < -5),:)));
intstuff.eparea = sum(ep_area);

forplotdvp = dVp_sumsub(T >= -5 & T <= 40,:);
forplotdvpm = dVpm_sumsub(T >= -5 & T <= 40);
forplotepm = mean_ep_b_eit(T >= -5 & T <= 40);
s = 0;
%for j = 1:size(dV_sigF_sumsub,2)
%    forfiltsumsub(1+s:50000+s,1) = detrend(dV_sigF_sumsub(:,j));
%    s = s + 50000;
%end

colorinc = 0;
figure;
set(gca,'FontSize',18)
hold on
for v = 1:size(EPall,1)
    plot(T,(detrend(EPall(v,:,12)))/1000,'color',[1-colorinc 0 0],'linewidth',3);
    colorinc = colorinc + 1/(size(EPall,1));
end
xlabel('Time (ms)','FontSize', 18);
ylabel('EP (mV)','FontSize', 18);
xlim(xlims);
title(sprintf('CAP over all trials on elec %d\n',12),'FontSize', 18);
hold off
drawnow;