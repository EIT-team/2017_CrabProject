HDR=ScouseTom_getHDR;
[~, fname]=fileparts(HDR.FileName);

Trigger=ScouseTom_TrigReadChn(HDR);
TT=ScouseTom_TrigProcess(Trigger,HDR);
Fs=HDR.SampleRate;
%% read the data
% get the voltages on the desired channels out of the EEG structure
Data = sread(HDR,inf,0) ;

Data(:,18:end) = [] ;

good_chn = [3 4 5 6 7 8 9 10 11 12] ; 

eit_inj_pairs = [8 9];
eit_freq = 225;

other_chn = 1:size(Data,2) ;

injtime = 360; %seconds
injnum = 1; %injections per second
j = 1;
k = 1;
g = 0;

%% Sort the Triggers

    for i = 1:length(TT.Stimulations{1})
        trigbois(j,k) = TT.Stimulations{1}(i);
        j = j+1;
        if j > injtime*injnum
            k = k + 1;
            j = 1;
        end
    end

xlims = [-5 40];

start_trial = 4 ; % trial we want to start with, for some reason the first few are fucked

T_trig = trigbois(1:injnum*injtime,1); % window in ms around event to view

tau_max = 250; % specify in ms

% find max timing between stims
Tmax = mean(floor((diff(T_trig)*1000)))/Fs;

tau = min([ tau_max Tmax]);

size_bin=floor(tau*Fs/1000); % convert to the number of samples this is equivalent to

%% Estimate the Carrier Frequency

% find first chunk of data
est_trig = (injtime*injnum)/2; % dont use the first one as sometimes its messed up
est_seg=detrend(Data(floor(T_trig(est_trig)*0.99):ceil(T_trig(est_trig+1)*1.01),:));

%find the carrier using this chunk
Fc_est = ScouseTom_data_GetCarrier(est_seg(:,eit_inj_pairs(2)),Fs);
Fc = round(Fc_est);

[bep aep] = butter(3,100/(Fs/2),'low');
[bepn aepn] = iirnotch(eit_freq/(Fs/2),(eit_freq/(Fs/2))/35);
DataF_EPlp = filtfilt(bep,aep,Data);
DataF_EP = filtfilt(bepn,aepn,DataF_EPlp);

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

%% spliting into each stim event
% make a Time vector - for plotting
T = (1:size_bin)*1000/Fs;
T = T - T(round(length(T)/2));

T_step = T(101)-T(100);

disp('Segmenting');
Data_seg=zeros(length(T_trig),size_bin,size(Data,2));
EPall=zeros(length(T_trig),size_bin,size(Data,2));
EITall=zeros(length(T_trig),size_bin,size(Data,2));

% loop through every stime event and take the data either side of stim
for iTrig=1:length(T_trig)
    Data_seg(iTrig,:,:)= Data((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:);
    EPall(iTrig,:,:)= DataF_EP((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:);
    EITall(iTrig,:,:)= DataF_EIT((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:);
end

%{
for n = 4:40

    clear hupper_test

    hupper_test = abs(hilbert(dV_sigF(:,n)));
    hupperf_test = filtfilt(FiltLPdV,hupper_test);
    figure
    plot(T,hupperf_test);
    hold on
    plot(T,dV_sigF(:,n));
    hold off
    xlim(xlims);
    drawnow
    
end
%}
%% Do simple coherent averaging 

% average all the EP chunks across repeats of EPs
EP_avg=detrend(squeeze(mean(EPall(start_trial:end,:,:),1)));
EIT_avg=detrend(squeeze(mean(EITall(start_trial:end,:,:),1)));

% find 
EIT_BV=squeeze(mean(EITall(start_trial:end, T > -100 & T < -30,:),2));
EIT_BVm=mean(EIT_BV);
EIT_BVstd=std(EIT_BV);

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
errorbar(EIT_BVm/1000,EIT_BVstd/1000);
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

%% EIT RandPhase

Y = Data_seg(4:end,:,eit_inj_pairs(1)-1)';

if mod(size(Y,2),2)==1
    Y = Y(:,2:end);
end

%% EIT from sum sub

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
%{
FiltLPdV = designfilt('lowpassfir', ...       % Response type
    'FilterOrder',500, ...            % Filter order
    'CutoffFrequency',40, ...    % Frequency constraints
    'DesignMethod','window', ...         % Design method
    'Window','blackmanharris', ...         % Design method options
    'SampleRate',Fs);               % Sample rate

hupperf = filtfilt(FiltLPdV,hupper);

figure;
plot(T,hupperf);
title('Low-pass Filtered Hilbert Transform Signal');
xlim(xlims);
drawnow;
%}
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

drawnow