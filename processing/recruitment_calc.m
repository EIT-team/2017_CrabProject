HDR=ScouseTom_getHDR;
[~, fname]=fileparts(HDR.FileName);

Trigger=ScouseTom_TrigReadChn(HDR);
TT=ScouseTom_TrigProcess(Trigger,HDR);
Fs=HDR.SampleRate;
%% read the data
% get the voltages on the desired channels out of the EEG structure
Data = sread(HDR,inf,0) ;

Data(:,18:end) = [] ;

good_chn = [3 4] ; 

injtime = 10; %seconds
injnum = 1; %injections per second

j = 1;
k = 1;
g = 0;

xlims = [-5 40];

%% Sort the Triggers

    for i = 1:length(TT.Stimulations{1})
        trigbois(i) = TT.Stimulations{1}(i);
    end
    
start_trial = 2 ; % trial we want to start with, for some reason the first few are fucked

%% Set this bad boy up

T_trig = trigbois(1:injnum*injtime); % window in ms around event to view

tau_max = 250; % specify in ms

% find max timing between stims
Tmax = mean(floor((diff(T_trig)*1000)))/Fs;

tau = min([ tau_max Tmax]);

size_bin=floor(tau*Fs/1000); % convert to the number of samples this is equivalent to

%% Estimate the Carrier Frequency

% find first chunk of data
est_trig = (injtime*injnum)/2; % dont use the first one as sometimes its messed up

N = 500;
FcEP = 100;

FiltEP = designfilt('lowpassfir', ...       % Response type
    'FilterOrder',N, ...            % Filter order
    'CutoffFrequency',FcEP, ...    % Frequency constraints
    'DesignMethod','window', ...         % Design method
    'Window','blackmanharris', ...         % Design method options
    'SampleRate',Fs);               % Sample rate

disp('Filtering for EPs');

DataF_EP = Data;%filtfilt(FiltEP,Data);
   
%% spliting into each stim event
% make a Time vector - for plotting
T = (1:size_bin)*1000/Fs;
T = T - T(length(T)/2);

T_step = T(101)-T(100);

% pre allocate the data

disp('Segmenting');
Data_seg=zeros(length(T_trig),size_bin,size(Data,2));
EPall=Data_seg;

% loop through every stime event and take the data either side of stim
for iTrig=1:length(T_trig)
    Data_seg(iTrig,:,:)= Data((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:);
    EPall(iTrig,:,:)= DataF_EP((T_trig(iTrig)-floor(size_bin/2):T_trig(iTrig)+ceil(size_bin/2)-1),:);
end



%% Do simple coherent averaging 

% average all the EP chunks across repeats of EPs
EP_avg=detrend(squeeze(mean(EPall(start_trial:end,:,:),1)));

EParea_t = zeros(size(T,2),3);
EParea_f = zeros(size(T,2),3);
EPamp_t = zeros(size(T,2),3);
EPamp_f = zeros(size(T,2),3);

%%
pstim_line_t(size(EP_avg(:,3)),1) = mean(EP_avg(12000:12400,3));

for l=1:5500 
    meow = l + 12000;
    EPdelta_t = EP_avg(meow,3)-mean(EP_avg(12000:12400,3));
    if meow >= 12680 && meow < 12980 && EPdelta_t <= 0
        EParea_t(meow,1) = T_step*EPdelta_t;
        EPamp_t(meow,1) = EP_avg(meow,3);
      
    end
    if meow >= 12980 && meow < 13342 && EPdelta_t <= 0
        EParea_t(meow,2) = T_step*EPdelta_t;
        EPamp_t(meow,2) = EP_avg(meow,3);
    end
    if meow >= 13342 && meow < 13829 && EPdelta_t <= 0
        EParea_t(meow,3) = T_step*EPdelta_t;
        EPamp_t(meow,3) = EP_avg(meow,3);
    end
end

pstim_line_f(size(EP_avg(:,4)),1) = mean(EP_avg(12000:12400,4));

for l=1:5500 
    meow = l + 12000;
    EPdelta_f = EP_avg(meow,4)-mean(EP_avg(12000:12400,4));
    if meow >= 12764 && meow < 13323 && EPdelta_f <= 0
        EParea_f(meow,1) = T_step*EPdelta_f;
        EPamp_f(meow,1) = EP_avg(meow,4);
    end
    if meow >= 13323 && meow < 13902 && EPdelta_f <= 0
        EParea_f(meow,2) = T_step*EPdelta_f;
        EPamp_f(meow,2) = EP_avg(meow,4);
    end
    if meow >= 13902 && meow < 14412 && EPdelta_f <= 0
        EParea_f(meow,3) = T_step*EPdelta_f;
        EPamp_f(meow,3) = EP_avg(meow,4);
    end
end

EPamp_min_t(1) = min(EPamp_t(:,1));
EPamp_min_t(2) = min(EPamp_t(:,2));
EPamp_min_t(3) = min(EPamp_t(:,3));

EPamp_min_f(1) = min(EPamp_f(:,1));
EPamp_min_f(2) = min(EPamp_f(:,2));
EPamp_min_f(3) = min(EPamp_f(:,3));

EParea_sum_t(1) = sum(EParea_t(:,1));
EParea_sum_t(2) = sum(EParea_t(:,2));
EParea_sum_t(3) = sum(EParea_t(:,3));

EParea_sum_f(1) = sum(EParea_f(:,1));
EParea_sum_f(2) = sum(EParea_f(:,2));
EParea_sum_f(3) = sum(EParea_f(:,3));

EPcv(1) = 4/(T(12764+find(EP_avg(12764:13323,4) <= min(EP_avg(12764:13323,4)))) - T(12680+find(EP_avg(12680:12980,3) <= min(EP_avg(12680:12980,3)))));
EPcv(2) = 4/(T(13323+find(EP_avg(13323:13902,4) <= min(EP_avg(13323:13902,4)))) - T(12980+find(EP_avg(12980:13342,3) <= min(EP_avg(12980:13342,3)))));
EPcv(3) = 4/(T(13902+find(EP_avg(13902:14412,4) <= min(EP_avg(13902:14412,4)))) - T(13342+find(EP_avg(13342:13829,3) <= min(EP_avg(13342:13829,3)))));

EPstruc(1).EP = Data_seg;
EPstruc(1).EP_avg = EP_avg;

figure;

subplot(2,1,1)
hold on
title(sprintf('Simple Coherent Averaging of EP on elec 3 \n'));
h1=plot(T,EP_avg(:,3),'linewidth',3);
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
title(sprintf('Simple Coherent Averaging of EP on elec 4\n'));
h1=plot(T,EP_avg(:,4),'linewidth',3);
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
