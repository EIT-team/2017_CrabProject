HDR=ScouseTom_getHDR;
[~, fname]=fileparts(HDR.FileName);

Trigger=ScouseTom_TrigReadChn(HDR);
TT=ScouseTom_TrigProcess(Trigger,HDR);
Fs=HDR.SampleRate;
%% read the data
% get the voltages on the desired channels out of the EEG structure
Data = sread(HDR,inf,0) ;

Data(:,18:end) = [] ;

good_chn = [3 4 5 6 8 9 10 11 12] ; 

eit_inj_pairs = [5 6 8 9 9 10 11 12 5 6];

other_chn = 1:size(Data,2) ;

injtime = 30; %seconds
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

start_trial = 2 ; % trial we want to start with, for some reason the first few are fucked

%% Set this bad boy up

for i = 1:length(eit_inj_pairs)/2

T_trig = trigbois(1:injnum*injtime,i); % window in ms around event to view

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
Fc_est = ScouseTom_data_GetCarrier(est_seg(:,inj_chn(1)),Fs);
Fc = round(Fc_est);

N = 500;
FcEP = 150;

FiltEP = designfilt('lowpassfir', ...       % Response type
    'FilterOrder',N, ...            % Filter order
    'CutoffFrequency',FcEP, ...    % Frequency constraints
    'DesignMethod','window', ...         % Design method
    'Window','blackmanharris', ...         % Design method options
    'SampleRate',Fs);               % Sample rate

disp('Filtering for EPs');

DataF_EP = Data;
BW = 125;
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

EPstruc(i).EP = Data_seg;
EPstruc(i).EP_avg = EP_avg;
EPstruc(i).BV = EIT_BV;

%%

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
% ylim([-4100 9000])

legend(h1,HDR.Label{good_chn})
%{
subplot(3,1,2);
title('Simple Coherent Averaging - EIT filtered and demod');
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

legend(h2,HDR.Label{good_chn})

%}
subplot(2,1,2)
errorbar(EIT_BVm/1000,EIT_BVstd/1000);
xlabel('Electrode');
ylabel('BV mV')
title('Boundary Voltages')

drawnow

%% EIT Sumsub
Y = EPall(4:end,:,eit_inj_pairs(i+g)-1)';

if mod(size(Y,2),2)==1
    Y = Y(:,2:end);
end

A = detrend(Y(:,1:2:end-2),'constant');
B = detrend(Y(:,2:2:end-2),'constant');
C = detrend(Y(:,3:2:end),'constant');


%% EP from sum sub

[b,a] = butter(5,200*2/Fs,'low');
EP = detrend(-(A+B)/2,'constant');
EP = filtfilt(b,a,EP);
EPm= mean(EP,2);

%% EIT from sum sub
dV_sig_orig =(A-2*B+C)/4; % kirills line ar fit way

dV_sigF=dV_sig_orig;

BW = 100;
Fc = 225;
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
    
dV_sigF=filtfilt(FiltBP,dV_sigF);

[dVdemod,hlower] = envelope(dV_sigF);

BV= mean(dVdemod(T > - 80 & T < -20,:));

dV=dVdemod-BV;
dVm=mean(dV,2);
dVp=100*(dV./BV);
dVpm=mean(dVp,2);

dVstruc(i).dV = dV;
dVstruc(i).dVp = dVp;

figure
subplot(2,1,1)
hold on
plot(T,(dV),'color',[0.7 0.7 0.7]);
plot(T,dVm,'linewidth',3);
title(sprintf('dV on elec %d\n',eit_inj_pairs(i+g)-1));
ylabel('uV');
xlabel('T ms');
hold off
xlim(xlims);
ylim([-50 50]);

subplot(2,1,2)
hold on
plot(T,(dVp),'color',[0.7 0.7 0.7]);
plot(T,dVpm,'linewidth',3);
title(sprintf('dVp on elec %d\n',eit_inj_pairs(i+g)-1));
ylabel('%');
xlabel('T ms');
hold off
xlim(xlims);
ylim([-0.2 0.2]);
drawnow

g = g + 1;

end

figure
subplot(3,2,1)
hold on
for i = 1:length(eit_inj_pairs)/2
    plot(T,detrend(squeeze(mean(EPstruc(i).EP(start_trial:end,:,3),1))),'linewidth',3);
end
title('EP on Electrode 3 over all trials')
ylabel('uV');
xlabel('T ms');
hold off
xlim([xlims]);
ylim([-14000 ep_y_max]);
legend({'Trial 1','Trial 2','Trial 3','Trial 4','Trial 5'});

subplot(3,2,2)
hold on
for i = 1:length(eit_inj_pairs)/2
    plot(T,detrend(squeeze(mean(EPstruc(i).EP(start_trial:end,:,4),1))),'linewidth',3);
end
title('EP on Electrode 4 over all trials')
ylabel('uV');
xlabel('T ms');
hold off
xlim([xlims]);
ylim([-14000 ep_y_max]);
legend({'Trial 1','Trial 2','Trial 3','Trial 4','Trial 5'});

subplot(3,2,3)
hold on
for i = 1:length(eit_inj_pairs)/2
    plot(T,detrend(squeeze(mean(EPstruc(i).EP(start_trial:end,:,5),1))),'linewidth',3);
end
title('EP on Electrode 5 over all trials')
ylabel('uV');
xlabel('T ms');
hold off
xlim([xlims]);
ylim([-14000 ep_y_max]);
legend({'Trial 1','Trial 2','Trial 3','Trial 4','Trial 5'});

subplot(3,2,4)
hold on
for i = 1:length(eit_inj_pairs)/2
    plot(T,detrend(squeeze(mean(EPstruc(i).EP(start_trial:end,:,6),1))),'linewidth',3);
end
title('EP on Electrode 6 over all trials')
ylabel('uV');
xlabel('T ms');
hold off
xlim([xlims]);
ylim([-14000 ep_y_max]);
legend({'Trial 1','Trial 2','Trial 3','Trial 4','Trial 5'});

subplot(3,2,5)
hold on
for i = 1:length(eit_inj_pairs)/2
    plot(T,detrend(squeeze(mean(EPstruc(i).EP(start_trial:end,:,12),1))),'linewidth',3);
end
title('EP on Electrode 12 over all trials')
ylabel('uV');
xlabel('T ms');
hold off
xlim([xlims]);
ylim([-14000 ep_y_max]);
legend({'Trial 1','Trial 2','Trial 3','Trial 4','Trial 5'});

subplot(3,2,6)
hold on
for i = 1:length(eit_inj_pairs)/2
    plot(T,detrend(squeeze(mean(EPstruc(i).EP(start_trial:end,:,13),1))),'linewidth',3);
end
title('EP on Electrode 13 over all trials')
ylabel('uV');
xlabel('T ms');
hold off
xlim([xlims]);
ylim([-14000 ep_y_max]);
legend({'Trial 1','Trial 2','Trial 3','Trial 4','Trial 5'});