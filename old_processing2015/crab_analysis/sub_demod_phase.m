function [T,mdV,sdV,mdPhi]=sub_demod_phase(Y,Fc,tstr)
%
    Fs = 16384;
    T = 1e3*((1:Fs/2)/Fs-.15);
    tClip = 700;
    thres = 100;

    if mod(size(Y,2),2)==1
        Y = Y(:,2:end);
    end

    A = detrend(Y(:,1:2:end-2),'constant');
    B = detrend(Y(:,2:2:end-2),'constant');
    C = detrend(Y(:,3:2:end),'constant');

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

    [b,a] = butter(5,200*2/Fs,'low');
    EP = detrend(-(A+B)/2,'constant');
    EP(T>-1.5 & T<1.5,:) = 0;
    EP = filtfilt(b,a,EP);
    
    [dV, dPhi] = hilbert_demod((A-2*B+C)/4,Fc,Fs);
    mu = mean(mean(dV));
    
     
    disp(sprintf('Mean BV = %.3f mV',mu/1e3));

%         dV = detrend(dV,'constant');
%         dV(tClip+1:end-tClip,:) = detrend(dV(tClip+1:end-tClip,:));
%         dV(T>-3 & T<3,:) = 0;
%         s = all(abs(dV(T>-10 & T<50,:))<thres);
% 
%         dPhi = detrend(dPhi,'constant');
%         dPhi(tClip+1:end-tClip,:) = detrend(dPhi(tClip+1:end-tClip,:));
%         dPhi(T>-1.5 & T<1.5,:) = 0;


    
    
    %    disp(sprintf('Removed %i noisy trials (threshold = %i uV)',nnz(~s),thres));

    %             figure('Position',[100,0,560,950],'PaperPositionMode','auto')
    %         
    %             subplot(3,1,1)
    %             plot(T,Y*1e-3)
    %             ylim([-70,70])
    %             xlim([-10,30])
    %             grid on
    %             xlabel('Time (ms)')
    %             ylabel('BV (mV)')  
    %             title(sprintf('mean BV = %.3f mV',mu*1e-3))
    %             
    %             subplot(3,1,2)
    %         % $$$     plot(T,EP,'color',.6*[1,1,1]), hold on
    %         % $$$ plot(T,mean(EP,2),'k','LineWidth',2), hold on    
    %             mEP = mean(EP,2);
    %             sEP = std(EP,0,2)/sqrt(size(EP,2));
    %             plot(T,mEP,'k','LineWidth',2), hold on
    %             plot(T,mEP+1.96*sEP,'--k')
    %             plot(T,mEP-1.96*sEP,'--k')
    %             ylim([-5e2,1e3])
    %             xlim([-10,30])
    %             grid on
    %             xlabel('Time (ms)')
    %             ylabel('EP (uV)')
    %             title(tstr,'interpreter','none')
    %         
    %             subplot(3,1,3)
    % $$$     plot(T,100*dV(:,s)/mu,'color',.6*[1,1,1]), hold on
    % $$$     plot(T,100*mean(dV(:,s),2)/mu,'k','LineWidth',2)
    mdV = mean(dV,2);
    mdPhi = mean(dPhi,2);
  	sdV = 100;%*std(dV(:,s),0,2)/(mu*sqrt(nnz(s)));
% $$$     hold on
% $$$     plot(T,100*dV(:,s)/mu,'Color',.6*[1,1,1])
    %             plot(T,mdV,'k','LineWidth',2), hold on
    %             plot(T,mdV+1.96*sdV,'--k')
    %             plot(T,mdV-1.96*sdV,'--k')
    %             xlim([-10,30])
    %             ylim([-.1,.1])
    %             grid on
    %             xlabel('Time (ms)')
    %             ylabel('dZ (%)')
