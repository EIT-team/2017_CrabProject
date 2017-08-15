function [T,mdV,sdV,mEP,mu]=sub_demod(Y,Fc,tstr,if_plot,hasTrig, thres, ppp)
%
    Fs = 16384;
    T = 1e3*((1:Fs/2)/Fs-.15);
    tClip = 700;
%     thres=50;
%     if (Fc>500) thres = 30;end;

    if mod(size(Y,2),2)==1
        Y = Y(:,1:end-1);
    end

    A = detrend(Y(:,1:2:end-2),'constant');
    B = detrend(Y(:,2:2:end-2),'constant');
    C = detrend(Y(:,3:2:end),'constant');
%      Y = Y(:,hasTrig(2:end-1));
%     phase = sum(Y.*repmat(Y(:,1),1,size(Y,2)));
%     pos = find(phase>0);
%     neg = find(phase<0);
%     pos = pos(1:min(length(pos),length(neg)));
%     neg = neg(1:min(length(pos),length(neg)));
%     
%     A = detrend(Y(:,pos(1:end-1)),'constant');
%     B = detrend(Y(:,neg(1:end-1)),'constant');
%     C = detrend(Y(:,pos(2:end)),'constant');
    
    k=size(A,2);
       
%     for i=1:k
%         A(:,i)=correct_cos_stim_by_same_fun(T,A(:,i),Fs,Fc,1);
%         B(:,i)=correct_cos_stim_by_same_fun(T,B(:,i),Fs,Fc,0);
%         C(:,i)=correct_cos_stim_by_same_fun(T,C(:,i),Fs,Fc,1);
%     end
    
     sig=(A-2*B+C)/4;
%     men=mean(sig((T>-0.1 & T<=2),:),2);
%     for i=1:size(sig,2)
%         sig((T>-0.1 & T<=2),i)=men;
%     end
%     for i=1:k
%         A(:,i)=correct_cos_stimulation(T,A(:,i),Fs,Fc);
%         B(:,i)=correct_cos_stimulation(T,B(:,i),Fs,Fc);
%         C(:,i)=correct_cos_stimulation(T,C(:,i),Fs,Fc);
%     end
    % outlier rejection
    rY = range(Y);
    ff_A = find(abs(range(A)-median(rY))>3*std(rY));
    ff_B = find(abs(range(B)-median(rY))>3*std(rY));
    ff_C = find(abs(range(C)-median(rY))>3*std(rY));
    ff = unique([ff_A,ff_B,ff_C]);

    if ~isempty(ff)
        disp(sprintf('Removed %i outliers out of %i',...
                     length(ff),size(A,2)));
        A(:,ff) = [];
        B(:,ff) = [];
        C(:,ff) = [];
    end

    [b,a] = butter(5,150*2/Fs,'low');
    EP = detrend(-(A+B)/2,'constant');
    EP(T>-1.5 & T<1.5,:) = 0;
    EP = filtfilt(b,a,EP);
    
   
    
    dV = hilbert_demod(sig,Fc,Fs);
   %dV = amdemod((A-2*B+C)/4,Fc,Fs);
    mu = mean(mean(dV));
    disp(sprintf('Mean BV = %.3f mV',mu/1e3));

    dV = detrend(dV,'constant');
    dV(tClip+1:end-tClip,:) = detrend(dV(tClip+1:end-tClip,:));
    dV(T>-ppp & T<ppp,:) = 0;
    s = all(abs(dV(T>-10 & T<50,:))<thres);

  
    
    disp(sprintf('Removed %i noisy trials (threshold = %i uV)',nnz(~s),thres));

    if (if_plot)
            figure('Position',[100,0,560,950],'PaperPositionMode','auto')
        
            subplot(3,1,1)
%             plot(T,A*1e-3); hold on;
%             plot(T,B*1e-3)
            plot(T,sig*1e-3);
            ylim([-70,70])
            xlim([-10,30])
            grid on
            xlabel('Time (ms)')
            ylabel('BV (mV)')  
            title(sprintf('mean BV = %.3f mV',mu*1e-3))
            
            subplot(3,1,2)
              plot(T,EP,'color',.6*[1,1,1]), hold on
          plot(T,mean(EP,2),'k','LineWidth',2), hold on    
end;
             mEP = mean(EP,2);
            sEP = std(EP,0,2)/sqrt(size(EP,2));
            
if (if_plot)
            plot(T,mEP,'k','LineWidth',2), hold on
            plot(T,mEP+1.96*sEP,'--k')
            plot(T,mEP-1.96*sEP,'--k')
            ylim([-5e2,1e3])
            xlim([-10,30])
            grid on
            xlabel('Time (ms)')
            ylabel('EP (uV)')
            title(tstr,'interpreter','none')
        
            subplot(3,1,3)
            plot(T,100*dV(:,s)/mu,'color',0.6*[1,1,1]), hold on
         plot(T,100*mean(dV(:,s),2)/mu,'k','LineWidth',2)
end;
    mdV = 100*mean(dV(:,s),2)/mu;
    sdV = 100*std(dV(:,s),0,2)/(mu*sqrt(nnz(s)));
% $$$     hold on
% $$$     plot(T,100*dV(:,s)/mu,'Color',.6*[1,1,1])
%             plot(T,mdV,'k','LineWidth',2), hold on
%             plot(T,mdV+1.96*sdV,'--k')
%             plot(T,mdV-1.96*sdV,'--k')
if (if_plot)
            xlim([-10,30])
            ylim([-.1,.1])
            grid on
            xlabel('Time (ms)')
            ylabel('dZ (%)')
end;