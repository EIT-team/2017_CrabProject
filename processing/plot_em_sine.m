HDR = ScouseTom_getHDR;
Fs = HDR.SampleRate;
Data = sread(HDR,Inf,0);

%%
range = [1:4];
inj = [1,2];

T = 1000*[1:length(Data)]/Fs;
[b_cap,a_cap]=butter(3,[70/(Fs/2) 140/(Fs/2)]);
%[b,a] = fir1(1000,[600/(Fs/2) 5100/(Fs/2)],'bandpass');
%freqz(b,a)
%fData = filtfilt(b_cap,a_cap,Data(:,range));
fData = filtfilt(b_cap,a_cap,Data(:,1)-Data(:,2));

[b_stim,a_stim]=butter(1,[5/(Fs/2), 15/(Fs/2)]);
Stimu = filtfilt(b_stim,a_stim,Data(:,range));
[pks_stim,locs_stim]=findpeaks(Stimu(:,1),'MinPeakProminence',5e4,'MinPeakDistance',8000); %500 if you have low stim

T_trig_idx = 1:length(locs_stim);
T_trig_idx= T_trig_idx(4:end-4);

figure
plot(Stimu(:,1))
hold on
plot(locs_stim,pks_stim,'o')
title('Peaks in Stim signal')

%% 

dz_pk_chn=2; % which channel to use in analysis of EIT injection

% high pass to remove stim and find carrier
[b_hp,a_hp]=butter(3,15/(Fs/2),'high');

Fc= ScouseTom_data_GetCarrier(filtfilt(b_hp,a_hp,Data(:,dz_pk_chn)),Fs);
% Fc=225;
BW=125;

[b_eit,a_eit]=butter(3,[(Fc-BW)/(Fs/2), (Fc+BW)/(Fs/2)]);
% [b,a]=butter(3,[725/(Fs/2), 1325/(Fs/2)]);
dZData = filtfilt(b_eit,a_eit,Data(:,range));
% dZData = dZData(:,inj(1))-dZData(:,inj(2));
[pks_dz,locs_dz]=findpeaks(dZData(:,dz_pk_chn),'MinPeakProminence',1e4,'MinPeakDistance',50);

T_trig_idx_CAP_Start=T_trig_idx(locs_stim(T_trig_idx) <locs_dz(10));
T_trig_idx_CAP_End=T_trig_idx(locs_stim(T_trig_idx) >locs_dz(end-10));
T_trig_idx_EIT= T_trig_idx( T_trig_idx > max(T_trig_idx_CAP_Start) & T_trig_idx < min(T_trig_idx_CAP_End));

figure
plot(dZData(:,dz_pk_chn))
hold on
plot(locs_dz,pks_dz,'o')
title('Peaks in EIT signal')

%%
figure
plot(Stimu(:,1))
hold on
h1=plot(locs_stim(T_trig_idx_CAP_Start),pks_stim(T_trig_idx_CAP_Start),'o');
h2=plot(locs_stim(T_trig_idx_EIT),pks_stim(T_trig_idx_EIT),'o');
h3=plot(locs_stim(T_trig_idx_CAP_End),pks_stim(T_trig_idx_CAP_End),'o');
legend([h1 h2 h3],'Start CAP','EIT','End CAP')
title('Peaks in Stim signal - used in averaging')
%%
T_trig=locs_stim(T_trig_idx_EIT);


%DO SUBTRACTION
dZData=dZData(:,1)-dZData(:,2);



dZData_demod = (abs(hilbert(dZData)));
BV = mean(dZData_demod(T_trig(1):T_trig(end),:));

[b_eit_hp,a_eit_hp]=butter(3,20/(Fs/2),'high');

dZData_demod=filtfilt(b_eit_hp,a_eit_hp,dZData_demod);

w = [1:Fs/5]-Fs/10;
Tw = 1000*w/Fs;

CAP=zeros(length(w),size(fData,2));
dZ = zeros(length(w),size(dZData_demod,2));
Stim=zeros(length(w),1);
figure;
for i = 1:floor(size(T_trig,1)/2)*2
    CAP = CAP + fData(T_trig(i)+w,:)/length(T_trig);
    dZ = dZ + detrend(dZData_demod(T_trig(i)+w,:)/length(T_trig));
    plot(detrend(dZData_demod(T_trig(i)+w,:)/length(T_trig)));hold on;
    Stim = Stim + Stimu(T_trig(i)+w,:)/length(T_trig);
end

%%
figure
subplot(4,1,1);
plot(Tw,Stim);
legend(num2str(range'))
ylabel('Stim,uV')
subplot(4,1,2);
plot(Tw,CAP);
ylabel('CAP,uV')
% xlabel('time,ms')
%xlim([-80,80])
grid on
subplot(4,1,3);
% plot(Tw,100*(dZ-dZ(find(Tw==-40),:))./BV,'linewidth',2);
plot(Tw,dZ)
ylabel('dZ,uV')
xlabel('time,ms')
subplot(4,1,4);
plot(Tw,100*(dZ-dZ(find(Tw==-40),:))./(BV),'linewidth',2);
% plot(Tw,dZ)
ylabel('dZ,%')
xlabel('time,ms')
% ylim([-3,3])
grid on
%%

% plot of average CAP only stils sines

% plot of CAPS only with mean without sines - start and finish compared

% as above but with EIT too


