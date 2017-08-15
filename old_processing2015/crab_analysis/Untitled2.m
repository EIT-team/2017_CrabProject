clear Y;
clear EP;
clear m;
EIT=29;
n=[21];
p=[3]
for j=4
    k=1;
    for i=1:length(n)
    [T,mdV,sdV,mEP,mu]=inj_elec_analysis(EIT,n(i),p(i),[15,6],[-2e3,7e3],[-0.1,0.1]);
    Y(j-3,k,:)=mdV;
    EP(j-3,k,:)=mEP;
    m(j-3,k,:)=[mu,mu];
    k=k+1;
    end
end

h.EP=EP;
h.BV=m;
h.dZ=Y;
h.EIT=EIT;
h.n=n;
save('channel5_6_dogs_dinner_hooks.mat','h');

clear k;
for i=1:size(h.BV,2)
    k(:,i,:)=(m(:,i,:));
   
end
clear G;
figure;
for j=1
    
 subplot(3,1,1)
 clear G;
        G(:,:)=EP(j,:,:)/1e3;
        plot(T,G,'LineWidth',2)
        xlim([-10,35])
        %ylim([-3,5])
        grid on
        xlabel('Time (ms)')
        ylabel('EP (mV)')
        title (['electrode #' num2str(j+3)])
        subplot(3,1,2);
        G(:,:)=Y(j,:,:);
        plot(T,G,'LineWidth',2)
        xlim([-10,35])
        ylim([-.2,.2])
        grid on
        xlabel('Time (ms)')
        ylabel('dZ (%)')
        title (['electrode #' num2str(j+3)])
        subplot(3,1,3);
        clear G; G(:,:)=k(:,:)/1e3;
        plot([0,1],G,'LineWidth',2)
        xlim([0,1])
        %ylim([-1,1])
        grid on
        xlabel('')
        ylabel('BV, mV')
        title (['electrode #' num2str(j+3)])
end