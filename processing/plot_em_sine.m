HDR = ScouseTom_getHDR;
Fs = HDR.SampleRate;
Data = sread(HDR,Inf,0);


%%
range = [1:2];
inj = [1,2];

T = 1000*[1:length(Data)]/Fs;
[b,a]=butter(3,[70/(Fs/2) 140/(Fs/2)]);
%[b,a] = fir1(1000,[600/(Fs/2) 5100/(Fs/2)],'bandpass');
%freqz(b,a)
fData = filtfilt(b,a,Data(:,range));

[b,a]=butter(1,[5/(Fs/2), 15/(Fs/2)]);
Stimu = filtfilt(b,a,Data(:,1));
[pks,locs]=findpeaks(Stimu,'MinPeakProminence',5e4,'MinPeakDistance',8000);

T_trig = locs(4:end-4);
w = [1:Fs/5]-Fs/10;
Tw = 1000*w/Fs;

figure
plot(Stimu)
hold on
scatter(locs,pks)

%[b,a]=butter(3,[100/(Fs/2), 350/(Fs/2)]);
[b,a]=butter(3,[725/(Fs/2), 1325/(Fs/2)]);
dZData = filtfilt(b,a,Data(:,range));
dZData = dZData(:,inj(1))-dZData(:,inj(2));
[pks,locs]=findpeaks(dZData,'MinPeakProminence',5e4,'MinPeakDistance',50);

T_trig(T_trig<locs(10))=[];
T_trig(T_trig>locs(end-10))=[];



figure
plot(dZData)
hold on
scatter(locs,pks)

dZData = (abs(hilbert(dZData)));

BV = mean(dZData(T_trig(1):T_trig(end)));


CAP=zeros(length(w),size(fData,2));
dZ = zeros(length(w),1);
for i = 1:floor(size(T_trig,1)/2)*2
    CAP = CAP + fData(T_trig(i)+w,:)/length(T_trig);
    dZ = dZ + detrend(dZData(T_trig(i)+w,:)/length(T_trig));
end

%%
figure
subplot(2,1,1);
plot(Tw,CAP);
ylabel('CAP,uV')
xlabel('time,ms')
%xlim([-80,80])
grid on
subplot(2,1,2);
plot(Tw,100*(dZ-dZ(find(Tw==-40)))/BV,'k','linewidth',2);
ylabel('dZ,%')
xlabel('time,ms')
ylim([-3,3])
grid on
%%
% 
% reg = [T>=0000&T<109300];
% figure;
% subplot(2,1,1)
% plot(T(reg),fData(reg,:));
% ylabel('CAP,uV')
% %xlim([9,9.3])
% subplot(2,1,2)
% 
% %legend('')
% 
% plot(T(reg),detrend(Stimu(reg)));
% xlabel('time,s')
% ylabel('Stim, uV')
% %xlim([9,9.3])
% 
% 
% %%
% 
% % figure;
% % subplot(2,1,1)
% % plot(T,fData);
% ylabel('CAP,uV')
% xlim([9,9.3])
% subplot(2,1,2)
% [b,a]=butter(1,[5/(Fs/2), 15/(Fs/2)]);
% legend('')
% Stimu = filtfilt(b,a,Data(:,1));
% plot(T(T>6&T<10),detrend(Stimu(T>6&T<10)));
% xlabel('time,s')
% ylabel('Stim, uV')
% xlim([9,9.3])
% 


