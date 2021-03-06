clear all
%% Read the data
HDR=ScouseTom_getHDR;
[~, fname]=fileparts(HDR.FileName);
Trigger=ScouseTom_TrigReadChn(HDR);
TT=ScouseTom_TrigProcess(Trigger,HDR); %Triggers
Fs=HDR.SampleRate; % Sampling Frequency
Data = sread(HDR,inf,0); % Read all data from all electrodes
Data(:,18:end) = []; % Truncate electrodes you are not using (speeds it up)
%% Settings - CHANGE THESE FOR YOUR SPECIFIC EXPERIMENTAL PROTOCOL
good_chn = [8 16]; % Channels for Analysis 
injtime = 20; % Total time of recording
injnum = 1; % Number of injections per second
start_trial = 2 ; % Typically we start on the second one to disregard anything weird with timing on the first
%% Sort the Stimulation Triggers
for i = 1:length(TT.Stimulations{1})
    trigbois(i) = TT.Stimulations{1}(i);
end
%% Figuring out Timing Windows
T_trig = trigbois(1:injnum*injtime); % window in ms around event to view
tau_max = 250; % specify in ms
Tmax = mean(floor((diff(T_trig)*1000)))/Fs; % find max timing between stims
tau = min([ tau_max Tmax]); % choose whichever is smallest
size_bin=floor(tau*Fs/1000); % convert to the number of samples this is equivalent to

%% Handing EP Data
[bepnf aepnf] = iirnotch(50/(Fs/2),(50/(Fs/2))/35); % 50 Hz notch filter
DataF_EP = filtfilt(bepnf,aepnf,Data); % apply 50 Hz notch filter

%% Segmenting Data
T = (1:size_bin)*1000/Fs; % make a time vector
T = T - T(round(length(T)/2));
T_step = T(101)-T(100); % calculate one step size of the time vector
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
%% Analyze the EPs
% Preallocate the area and amplitude data arrays
EParea_t = zeros(size(T,2),3);
EParea_f = zeros(size(T,2),3);
EPamp_t = zeros(size(T,2),3);
EPamp_f = zeros(size(T,2),3);
% Calculate the average value of the data before stim event for the first
% channel
pstim_line_t(size(EP_avg(:,good_chn(1))),1) = mean(EP_avg(T > - 5 & T < -1,good_chn(1)));
% Calculate the areas of the different sections
for l=1:5500 % Only calculating from 0 to 50 milliseconds after stim -- YOU WILL NEED TO CHANGE IF YOU NEED BIGGER VIEW 
    i = l + 12000; % Dummy variable for keeping index right
    EPdelta_t = EP_avg(i,good_chn(1))-mean(EP_avg(12000:12400,good_chn(1))); % Calculate the amplitude of the change
    % This splits into the different areas based on time. It could be
    % possible to calculate these time points, but we have just empirically
    % determined them. 
    if i >= 12680 && i < 12980 && EPdelta_t <= 0 
        EParea_t(i,1) = T_step*EPdelta_t;
        EPamp_t(i,1) = EP_avg(i,good_chn(1));
    end
    if i >= 12980 && i < 13342 && EPdelta_t <= 0
        EParea_t(i,2) = T_step*EPdelta_t;
        EPamp_t(i,2) = EP_avg(i,good_chn(1));
    end
    if i >= 13342 && i < 13829 && EPdelta_t <= 0
        EParea_t(i,3) = T_step*EPdelta_t;
        EPamp_t(i,3) = EP_avg(i,good_chn(1));
    end
end
% All this is the same as above just for the second channel of interest. If
% you have more channels of interest you will need to copy/paste this
% section and create more variable names, etc.
pstim_line_f(size(EP_avg(:,good_chn(2))),1) = mean(EP_avg(12000:12400,good_chn(2)));
for l=1:5500 
    i = l + 12000;
    EPdelta_f = EP_avg(i,good_chn(2))-mean(EP_avg(12000:12400,good_chn(2)));
    if i >= 12764 && i < 13323 && EPdelta_f <= 0
        EParea_f(i,1) = T_step*EPdelta_f;
        EPamp_f(i,1) = EP_avg(i,good_chn(2));
    end
    if i >= 13323 && i < 13902 && EPdelta_f <= 0
        EParea_f(i,2) = T_step*EPdelta_f;
        EPamp_f(i,2) = EP_avg(i,good_chn(2));
    end
    if i >= 13902 && i < 14412 && EPdelta_f <= 0
        EParea_f(i,3) = T_step*EPdelta_f;
        EPamp_f(i,3) = EP_avg(i,good_chn(2));
    end
end
% Finding the minimum amplitudes of all the different sections for first
% channel
EPamp_min_t(1) = min(EPamp_t(:,1));
EPamp_min_t(2) = min(EPamp_t(:,2));
EPamp_min_t(3) = min(EPamp_t(:,3));
% Finding the minimum amplitudes of all the different sections for second
% channel
EPamp_min_f(1) = min(EPamp_f(:,1));
EPamp_min_f(2) = min(EPamp_f(:,2));
EPamp_min_f(3) = min(EPamp_f(:,3));
% Finding the total area for each different section for the first channel
EParea_sum_t(1) = sum(EParea_t(:,1));
EParea_sum_t(2) = sum(EParea_t(:,2));
EParea_sum_t(3) = sum(EParea_t(:,3));
% Finding the total area for each different section for the second channel
EParea_sum_f(1) = sum(EParea_f(:,1));
EParea_sum_f(2) = sum(EParea_f(:,2));
EParea_sum_f(3) = sum(EParea_f(:,3));
% Finding the conduction velocity for each different section for the first
% channel -- NOTE: The 4 in the beginning is because we have 4 mm distance
% between electrodes, you may need to change
EPcv_t(1) = (4*(good_chn(1)-1))/(T(12764+find(EP_avg(12764:13323,good_chn(1)) <= min(EP_avg(12764:13323,good_chn(1))))));
EPcv_t(2) = (4*(good_chn(1)-1))/(T(13323+find(EP_avg(13323:13902,good_chn(1)) <= min(EP_avg(13323:13902,good_chn(1))))));
EPcv_t(3) = (4*(good_chn(1)-1))/(T(13902+find(EP_avg(13902:14412,good_chn(1)) <= min(EP_avg(13902:14412,good_chn(1))))));
% Finding the conduction velocity for each different section for the second
% channel -- NOTE: The 4 in the beginning is because we have 4 mm distance
% between electrodes, you may need to change
EPcv_f(1) = (4*(good_chn(2)-1))/(T(12764+find(EP_avg(12764:13323,good_chn(2)) <= min(EP_avg(12764:13323,good_chn(2))))));
EPcv_f(2) = (4*(good_chn(2)-1))/(T(13323+find(EP_avg(13323:13902,good_chn(2)) <= min(EP_avg(13323:13902,good_chn(2))))));
EPcv_f(3) = (4*(good_chn(2)-1))/(T(13902+find(EP_avg(13902:14412,good_chn(2)) <= min(EP_avg(13902:14412,good_chn(2))))));
% Store data if you need it for later
EPstruc(1).EP = Data_seg;
EPstruc(1).EP_avg = EP_avg;
%{
% Frequency Analysis of the CAP -- NOTE: Y is the raw data from the CAP
L = 25000; % Length of the data array
fft_cap = fft(Y); % Take fft
P2 = abs(fft_cap/L); % Only want the power
P1 = P2(1:(L/2)+1); % Reshape
P1(2:end-1) = 2*P1(2:end-1); % Rescale
f = Fs * (0:(L/2))/L; % Find frequency axis based on sampling frequency
figure
plot(f,P1);
%}
%% Plotting Everything
xlims = [-5 40]; % Set x-limits
figure;
subplot(2,1,1)
hold on
title(sprintf('Simple Coherent Averaging of EP on first channel\n'));
h1=plot(T,EP_avg(:,good_chn(1)),'linewidth',3);
plot(T,pstim_line_t,'--','linewidth',2);
bar(T(12001:17500),EParea_t(12001:17500,1)/T_step);
bar(T(12001:17500),EParea_t(12001:17500,2)/T_step);
bar(T(12001:17500),EParea_t(12001:17500,3)/T_step);
xlabel('Time ms');
ylabel('EP uv');
xlim(xlims)
ep_y_max = round(max(max((1000+EP_avg( T > 2 & T <xlims(2),good_chn)))),-2);
ep_y_min = -round(max(max((1000-EP_avg( T > 2 & T <xlims(2),good_chn)))),-2);
ylim([ep_y_min ep_y_max])

subplot(2,1,2)
hold on
title(sprintf('Simple Coherent Averaging of EP on second channel\n'));
h1=plot(T,EP_avg(:,good_chn(2)),'linewidth',3);
plot(T,pstim_line_f,'--','linewidth',2);
bar(T(12001:17500),EParea_f(12001:17500,1)/T_step);
bar(T(12001:17500),EParea_f(12001:17500,2)/T_step);
bar(T(12001:17500),EParea_f(12001:17500,3)/T_step);
xlabel('Time ms');
ylabel('EP uv');
xlim(xlims)
ep_y_max = round(max(max((1000+EP_avg( T > 2 & T <xlims(2),good_chn)))),-2);
ep_y_min = -round(max(max((1000-EP_avg( T > 2 & T <xlims(2),good_chn)))),-2);
ylim([ep_y_min ep_y_max])

forplotEP = EP_avg(T >= -5 & T <= 40,good_chn(1));
forplotT = T(T >= -5 & T <= 40)';
