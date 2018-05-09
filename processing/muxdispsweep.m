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

eit_cur = 10; %eit current in uA

injtime = 60; %seconds
injnum = 5; %injections per second
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
Fc_est = ScouseTom_data_GetCarrier(est_seg(:,eit_inj_pairs(i+g)),Fs);
Fc = round(Fc_est);

disp('Filtering for EPs');

[bep aep] = butter(3,100/(Fs/2),'low');
[bepn aepn] = iirnotch(eit_freq/(Fs/2),(eit_freq/(Fs/2))/35);
DataF_EPlp = filtfilt(bep,aep,Data);
DataF_EP = filtfilt(bepn,aepn,DataF_EPlp);

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
T = T - T(round(length(T)/2));

T_step = T(101)-T(100);

% pre allocate the data
%{
Data_clean = zeros(size(Data,1),1);

    p = 2;
    
    trial_max = max(Data(((T_trig(10)-2800):T_trig(10)-800),7));
    trial_min = min(Data(((T_trig(10)-2800):T_trig(10)-800),7));
    
    for w = 1:size(Data,1)-1
       
        delta_V(p) = Data(w+1,7) - Data(w,7);
        cur_V(p) = Data(w,7);
        
        if delta_V(p) > 2000 || cur_V(p) > trial_max || cur_V(p) < trial_min
            Data_clean(w) = cur_V(p-1); 
        end
        if cur_V(p) <= trial_max && cur_V(p) >= trial_min
            Data_clean(w) = cur_V(p);
            p = p + 1;
        end
    end
%}
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
EPstruc(i).BVrecmean = mean(EIT_BV(:,eit_inj_pairs(i+g)-1));
EPstruc(i).BVlasteitmean = mean(EIT_BV(:,eit_inj_pairs(i+g+1)));
EPstruc(i).EPmaxf = min(EP_avg(12600:20000,eit_inj_pairs(i+g)-1));
EPstruc(i).EPmaxeit = min(EP_avg(12600:20000,eit_inj_pairs(i+g+1)));
EPstruc(i).EPcv = T(12600+find(EP_avg(12600:20000,4) <= min(EP_avg(12600:20000,4)))) - T(12600+find(EP_avg(12600:20000,3) <= min(EP_avg(12600:20000,3))));
t_cent_indx = find(EP_avg(12600:20000,eit_inj_pairs(i+g+1)) <= EPstruc(i).EPmaxeit);
EPstruc(i).EPeiteltime = T(12600+t_cent_indx);
EPstruc(i).nervez = (mean(EIT_BV(:,eit_inj_pairs(i+g)))-mean(EIT_BV(:,eit_inj_pairs(i+g+1))))/(eit_cur);
EPstruc(i).backcur_fivefour = (mean(EIT_BV(:,eit_inj_pairs(i+g)-1))-mean(EIT_BV(:,eit_inj_pairs(i+g))))/(EPstruc(i).nervez);
EPstruc(i).backcur_sixseven = (mean(EIT_BV(:,eit_inj_pairs(i+g+1)))-mean(EIT_BV(:,eit_inj_pairs(i+g+1)+1)))/(EPstruc(i).nervez);

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

%% EIT Sumsub
Y = Data_seg(4:end,:,eit_inj_pairs(i+g)-1)';

if mod(size(Y,2),2)==1
    Y = Y(:,2:end);
end

A = detrend(Y(:,1:2:end-2),'constant');
B = detrend(Y(:,2:2:end-2),'constant');
C = detrend(Y(:,3:2:end),'constant');

%% EIT from sum sub

dV_sig_orig = (A-2*B+C)/4; % kirills line ar fit way

BW = 100;
Fc = eit_freq;
N=1000;
F6dB1=Fc-BW;
F6dB2=Fc+BW;

[bbp,abp] = butter(3,[F6dB1 F6dB2]/(Fs/2));

dV_sigF=filtfilt(bbp,abp,dV_sig_orig);

hupper = abs(hilbert(dV_sigF));

FiltLPdV = designfilt('lowpassfir', ...       % Response type
    'FilterOrder',500, ...            % Filter order
    'CutoffFrequency',40, ...    % Frequency constraints
    'DesignMethod','window', ...         % Design method
    'Window','blackmanharris', ...         % Design method options
    'SampleRate',Fs);               % Sample rate

hupperf = filtfilt(FiltLPdV,hupper);

BV_u= mean(hupperf(T > - 80 & T < -20,:));

dV=hupperf-BV_u;
dVm=mean(dV,2);

dVp=100*(dV./BV_u);
dVpm=mean(dVp,2);

pstim_line(size(dVm),1) = mean(dVm(12000:12400));

dVa = [];

for l=1:5500
    dVdelta = dVm(l+12000)-mean(dVm(12000:12400));
    dVa(l) = T_step*dVdelta;
end


dVstruc(i).dV = dV;
dVstruc(i).dVp = dVp;
dVstruc(i).dVa = sum(dVa);
dVstruc(i).dVpmmin = min(dVpm(12600+t_cent_indx-200:12600+t_cent_indx+200)); %finding max dZ within 2 ms window from CAP on last EIT electrode
dVstruc(i).dVpmintime = T(12600+t_cent_indx-200+find(dVpm(12600+t_cent_indx-200:12600+t_cent_indx+200) <= dVstruc(i).dVpmmin));
dVstruc(i).dVmmin = min(dVm(12600+t_cent_indx-200:12600+t_cent_indx+200)); %finding max dZ within 2 ms window from CAP on last EIT electrode
dVstruc(i).dVmintime = T(12600+t_cent_indx-200+find(dVm(12600+t_cent_indx-200:12600+t_cent_indx+200) <= dVstruc(i).dVmmin));

figure
subplot(3,1,1)
hold on
plot(T,(dV),'color',[0.7 0.7 0.7]);
plot(T,dVm,'linewidth',3);
plot(T,pstim_line,'--','linewidth',2);
bar(T(12001:17500),dVa/T_step);
title(sprintf('dV on elec %d\n',eit_inj_pairs(i+g)-1));
ylabel('uV');
xlabel('T ms');
hold off
xlim(xlims);

subplot(3,1,2)
hold on
plot(T,(dVp),'color',[0.7 0.7 0.7]);
plot(T,dVpm,'linewidth',3);
title(sprintf('dVp on elec %d\n',eit_inj_pairs(i+g)-1));
ylabel('%');
xlabel('T ms');
hold off
xlim(xlims);

subplot(3,1,3)
hold on
bar(T(12001:17500),dVa);
title(sprintf('Area of dV on elec %d\n',eit_inj_pairs(i+g)-1));
ylabel('uV*ms');
xlabel('T ms');
hold off
xlim(xlims);

drawnow

g = g + 1;

end