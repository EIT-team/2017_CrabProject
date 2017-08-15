clear X;
% figure;
% 
% start=3;
% for i=4:8
% %     for j=start:11;
% %         if (j==i && i~=7) ch=[i,i+1];
% %         else if (j==i && i==7) ch=[i,i+2];
% %             else ch=[j,12];
% %             end
% %         end
%         ch=[i,i+1]
%         [T,mdV,sdV]=inj_elec_analysis(29,4,i-3,ch,[-2e3,7e3],[-0.1,0.1]);
%         X(i-3,:)=mdV;
%         subplot(5,1,i-3);
%         plot(T,mdV,'LineWidth',2), hold on
%         plot(T,mdV+1.96*sdV,'--k')
%         plot(T,mdV-1.96*sdV,'--k')
%         xlim([0,25])
%         ylim([-.04,.04])
%         grid on
%         xlabel('Time (ms)')
%         ylabel('dZ (%)')
%         title(['I ' num2str(i) '-' num2str(i+1) ', channel ' num2str(ch(1)) '-' num2str(ch(2))])
% %     end
% end
% 

for p=1:3

    clear Re;
    clear Im;
    clear dV;
    clear Base;
    j=0;

    if (p==1)
        ref=[5,15];
        chan=[4,7];
    elseif p==2
        ref=[6,16];
        chan=[4,8];
    elseif p==3
        ref=[7,17];
        chan=[4,9];
    end
   % chan=[4,19];

  %  for i=[7,8,9,11,21,22,24] % 10uA
   
  %  for i=[14,18,19] % 5uA
   for i=[10,12,13,15,16,17,20,23] % 20uA
   
    j=j+1;
    oo=10;
    %if(i>=11) oo=19;end;
    [T,mdV1,sdV,dphi_curr,Curr]=inj_elec_analysis_phase(30,i,p,ref,[-2e3,7e3],[-0.1,0.1]);
    
    [T,mdV,sdV,dphi_volt,Volt]=inj_elec_analysis_phase(30,i,p,chan,[-2e3,7e3],[-0.1,0.1]);  
   
  %  [T,mdV,sdV,dphi_volt,Volt]=inj_elec_analysis_phase(32,i,p,[4,7],[-2e3,7e3],[-0.1,0.1]);
    
    Volt = detrend(Volt,'constant');
    k=1;
    for l=1:size(Volt,2)
         Ba(k)=phase_diff( Volt(:,l),Curr(:,l));
         k=k+1;
    end
    Base=mean(Ba);
    BV=mean(max(Volt));
    dphi=dphi_volt-dphi_curr;
    dphi=detrend(dphi,'constant');
    dphi(700+1:end-700,:) = detrend(dphi(700+1:end-700,:));
    %dphi=detrend(dphi,'constant');
    Phi=dphi+Base;
    Re(j,:)=mdV.*cos(Phi);
    Im(j,:)=mdV.*sin(Phi);
    dV(j,:)=mdV;
    b_Re=mean(Re(j,:));
    b_Im=mean(Im(j,:));
    b_dV=mean(dV(j,:));
    dV(j,:)=100*((dV(j,:)-b_dV)./b_dV);
    dV(j,:)=dV(j,:)-mean(dV(j,2000:4000));
    dV(j,T<2 & T>-1)=0;
    
    Re(j,:)=100*((Re(j,:)-b_Re)./b_Re);
    Re(j,:)=Re(j,:)-mean(Re(j,2000:4000));
    Re(j,T<2 & T>-1)=0;
    
    Im(j,:)=100*((Im(j,:)-b_Im)./b_Im);
    Im(j,:)=Im(j,:)-mean(Im(j,2000:4000));
    Im(j,T<2 & T>-1)=0;
   
end
%     dA_p= 100*(1-(cos(base)./cos(dPhi)));
%     dA_p=detrend(dA_p,'constant');
%     dA_p=dA_p-dA_p(1000);

h.dV=dV;
h.Re=Re;
h.Im=Im;
h.T=T;

save(['Control_artifact_amp_225Hz_20uA_p' num2str(p) '.mat'],'h');

figure;
subplot(3,1,1);
    plot(T,dV(1:end,:),'Color',0.6*[1,1,1],'LineWidth',1)
    hold on
    plot(T,mean(dV(1:end,:)),'r','LineWidth',2)
    xlim([-10,60])
    ylim([-.15,.15])
    grid on;
    xlabel('Time (ms)')
    ylabel('d|Z|,%')
    subplot(3,1,2);
    plot(T,Re,'Color',0.6*[1,1,1],'LineWidth',1)
    hold on
    plot(T,mean(Re([1:end],:)),'r','LineWidth',2)
    xlim([-10,60])
    ylim([-0.15,.15])
    grid on;
    xlabel('Time (ms)')
    ylabel('dR,%')
    subplot(3,1,3);
    plot(T,Im,'Color',0.6*[1,1,1],'LineWidth',1)
    hold on
    plot(T,mean(Im),'r','LineWidth',2)
    xlim([-10,60])
    ylim([-5,5])
    grid on;
    xlabel('Time (ms)')
    ylabel('dX,%')
end
