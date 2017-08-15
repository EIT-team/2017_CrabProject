function Process_multichanel(EIT,p,n,num,title)


for j=4:(num+3)
    k=1;
    for i=1:length(n)
    %if(n(i)==3 || n(i)==21 || n(i)==20) p=2; end;
    l=j;
    ref=12;
    if (ismember(n(i),[11:24]))
        ref=19;
    end
    if (ismember(n(i),[7:24]))
        if(j==5) 
            l=15; end;
        if(j==6) 
            l=16; end;
        if(j==7) 
            l=17; end;
    end
    [T,mdV,sdV,mEP,mu]=inj_elec_analysis(EIT,n(i),p,[l,ref],[-2e3,7e3],[-0.1,0.1]);
    Y(j-3,k,:)=mdV;
    EP(j-3,k,:)=mEP;
    m(j-3,k,:)=[mu,mu];
    k=k+1;
    end
end
h.dZ=Y;
h.EP=EP;
h.BV=m;
h.EIT=29;
h.n=n;
%save ('225_HOOKS_ref12','h');

%load '225_HOOKS_ref12.mat'


Y=h.dZ;EP=h.EP;m=h.BV;
clear k;
for i=1:size(h.BV,2)
    k(:,i,:)=(m(:,i,:))./m(p+1,i,1);
end
clear G;

figure;
for j=1:num
    if (j>p+1) k(j,:,:)=-k(j,:,:); end;
 subplot(3,num,j)
 clear G;
        G(:,:)=EP(j,:,:)/1e3;
        plot(T,G,'LineWidth',2)
        xlim([-10,35])
        ylim([-3,10])
        grid on
        xlabel('Time (ms)')
        ylabel('EP (mV)')
        title (['electrode #' num2str(j+3)])
        subplot(3,num,j+num);
        clear G; G(:,:)=Y(j,:,:);
        plot(T,G,'LineWidth',2)
        xlim([-10,35])
        ylim([-.1,.1])
        grid on
        xlabel('Time (ms)')
        ylabel('dZ (%)')
        title (['electrode #' num2str(j+3)])
        subplot(3,num,j+num*2);
        clear G; G(:,:)=k(j,:,:);
        plot([0,1],G','LineWidth',2)
        xlim([0,1])
        ylim([-1,1])
        grid on
        xlabel('')
        ylabel('BV/BV_inj')
        title (['electrode #' num2str(j+3)])
end

