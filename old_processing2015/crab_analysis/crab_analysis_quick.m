function crab_analysis(Nerve_number, EIT_number, Ch_number, Fc)thres = 100;all_ch = [1:19,32];Fs = 16384;info.Fs = Fs;info.Stim_Rate = 2;info.Trig_Chan = length(all_ch);T = 1e3*((1:Fs/2)/Fs-.15);tClip = 700;if ~exist('mat','dir')    mkdir('mat')endfname = sprintf('mat/N%03i_EIT%03i.mat',Nerve_number,EIT_number);if exist(fname,'file')    load(fname,'Data','TrigData')else        EEG = pop_biosig(sprintf('N%03i_EIT%03i.bdf',Nerve_number,EIT_number));    EEG.data = EEG.data(all_ch,:);    EEG.data(end,:) = [0,diff(EEG.data(end,:))];    Data = fix_phase_jump(EEG.data',info);    TrigData = reshape(Data(:,end),Fs/2,[]);    save(fname,'Data','TrigData','all_ch')endhasTrig = find(any(abs(TrigData(T>-2 & T<2,:))>1e3));disp(sprintf('Found %i triggers',length(hasTrig)))if any(diff(hasTrig)>1)    warning('*******TRIGGER SKIPPED*********')enddisp(num2str(std(Data(:,1:end-1))))% figure('Position',get(0,'ScreenSize'))% subplot(1,3,1)% plot(T,TrigData(:,hasTrig))% xlim([-5,5])% subplot(1,3,2)% Y = reshape(Data(:,4),Fs/2,[]);% imagesc(detrend(Y(:,hasTrig(2:end-1)),'constant'))% subplot(1,3,3)% pwelch(Data(:,4),Fs,.95,[],Fs)% xlim([0,2e3])% % % for iCh = 1:length(Ch_number)%     tstr = sprintf('Nerve %03i, EIT %03i, Channel %i (Fc = %iHz)',...%                    Nerve_number,EIT_number,Ch_number(iCh),Fc);%     disp(tstr)%     %     sel_ch = find(all_ch==Ch_number(iCh));%     Y = reshape(Data(:,sel_ch),Fs/2,[]);    %     Y = Y(:,hasTrig(2:end-1));%     phase = sum(Y.*repmat(Y(:,1),1,size(Y,2)));%     pos = find(phase>0);%     neg = find(phase<0);%     pos = pos(1:min(length(pos),length(neg)));%     neg = neg(1:min(length(pos),length(neg)));%     %     A = detrend(Y(:,pos(1:end-1)),'constant');%     B = detrend(Y(:,neg(1:end-1)),'constant');%     C = detrend(Y(:,pos(2:end)),'constant');% %     % outlier rejection%     rY = range(Y);%     ff_A = find(abs(range(A)-median(rY))>3*std(rY));%     ff_B = find(abs(range(B)-median(rY))>3*std(rY));%     ff_C = find(abs(range(C)-median(rY))>3*std(rY));%     ff = unique([ff_A,ff_B,ff_C]);%     if ~isempty(ff)%         disp(sprintf('Removed %i outliers out of %i',...%                      length(ff),size(A,2)))%         A(:,ff) = [];%         B(:,ff) = [];%         C(:,ff) = [];%     end% %     [b,a] = butter(5,200*2/Fs,'low');%     EP = detrend(-(A+B)/2,'constant');%     EP(T>-1.5 & T<1.5,:) = 0;%     EP = filtfilt(b,a,EP);% %     dV = hilbert_demod(A-2*B+C,Fc,Fs);%     mu = mean(mean(dV));%     dV = detrend(dV,'constant');%     dV(tClip+1:end-tClip,:) = detrend(dV(tClip+1:end-tClip,:));%     dV(T>-1.5 & T<1.5,:) = 0;%     s = all(abs(dV(T>-10 & T<50,:))<thres);% %     disp(sprintf('Removed %i noisy trials (threshold = %i uV)',nnz(~s),thres))% %     figure('Position',[100,0,560,950],'PaperPositionMode','auto')%     subplot(2,1,1)%     plot(T,EP,'color',.6*[1,1,1]), hold on%     plot(T,mean(EP,2),'k','LineWidth',2)%     ylim([-1e3,1e3])%     xlim([-10,30])%     grid on%     xlabel('Time (ms)')%     ylabel('EP (uV)')% %     subplot(2,1,2)%     plot(T,100*dV(:,s)/mu,'color',.6*[1,1,1]), hold on%     plot(T,100*mean(dV(:,s),2)/mu,'k','LineWidth',2)%     xlim([-10,30])%     ylim([-.1,.1])%     grid on%     xlabel('Time (ms)')%     ylabel('dZ (%)')%     title(tstr)%     saveas(gcf,sprintf('N%03i_EIT%03i_Ch%i.png',Nerve_number,EIT_number,Ch_number(iCh)))% end