function [Tp,mdV,sdV,mEP,mu,Fc]=inj_elec_analysis(Nerve_number, EIT_number, Pair_number, Inj, y1, y2, if_plot, thres, ppp)
%

if nargin==4
    y1 = [-5e2,2e3];
    y2 = [-.1,.1];
end

load( sprintf('D:/CRAB/Summary_all/mat/N%03i_EIT%03i.mat',Nerve_number,EIT_number) );

Fs = 16384;
T = 1e3*((1:Fs/2)/Fs-.15);

thress=1e4;
if (Nerve_number==30) thress=1e3;end;

hasTrig = find(any(abs(TrigData(T>-2 & T<2,:))>thress));
% $$$ disp(sprintf('Found %i triggers',length(hasTrig)))
% $$$ if length(unique(diff(hasTrig)))>1
% $$$     warning('********* Trigger not detected in at least 1 trial ***********')
% $$$ end

pairChanges = [0,find(diff(hasTrig)>1),length(hasTrig)-1];
selPair = false(size(TrigData));
selPair(:,hasTrig(pairChanges(Pair_number)+2):...
          hasTrig(pairChanges(Pair_number+1)-1)) = true;
%Data([pairChanges,pairChanges+1],:) =[]; 
Data = Data(selPair(:),:);
disp(num2str([pairChanges(Pair_number)+2,pairChanges(Pair_number+1)-1]))

if length(Inj)==1
    Xinj = Data(:,all_ch==Inj(1));
    [Pxx,f] = pwelch(Xinj(:),Fs,.95,[],Fs);
    fsel = find(f>125);
    [Pmax,ind] = max(Pxx(fsel));
    Fc = f(fsel(ind));
    disp(sprintf('Detected: Fc = %i Hz',Fc))
    tstr = sprintf('[%i Hz] N%03i_EIT%03i_Pair%i: V%i',...
                   Fc,Nerve_number,EIT_number,Pair_number,Inj(1));
    fname = sprintf('N%03i_%04iHz_EIT%03i_Pair%i_V%i.png',...
                    Nerve_number,Fc,EIT_number,Pair_number,Inj(1));
else
    Xinj = Data(:,all_ch==Inj(1))-Data(:,all_ch==Inj(2));
    [Pxx,f] = pwelch(Xinj(:),Fs,.95,[],Fs);
    fsel = find(f>200);
    [Pmax,ind] = max(Pxx(fsel));
    Fc = f(fsel(ind));
    disp(sprintf('Detected: Fc = %i Hz',Fc))
    tstr = sprintf('[%i Hz] N%03i_EIT%03i_Pair%i: V%i-V%i',...
                   Fc,Nerve_number,EIT_number,Pair_number,Inj(1),Inj(2));
    fname = sprintf('N%03i_%04iHz_EIT%03i_Pair%i_V%i-V%i.png',...
                   Nerve_number,Fc,EIT_number,Pair_number,Inj(1),Inj(2));
end
Xinj = reshape(Xinj,Fs/2,[]);


% $$$ figure
% $$$ subplot(1,2,1)
% $$$ plot([std(X1)',std(X2)',std(X3)',std(X4)'])
% $$$ subplot(1,2,2)
% $$$ plot([std(X2)',std(X3)',std(X2-X3)'])
% $$$ waitforbuttonpress

k = size(Xinj,2);

% $$$ b23 = [ones(k,1),(0:k-1)']\std(Xinj1-Xinj2)';
% $$$ b2 = [ones(k,1),(0:k-1)']\std(Xinj1)';
% $$$ b3 = [ones(k,1),(0:k-1)']\std(Xinj2)';
% $$$ disp(sprintf('Slope: %.2f %.2f %.2f',b2(2),b3(2),b23(2)))

phase = sum(Xinj.*repmat(Xinj(:,1),1,k));
if length(unique([diff(find(phase>0)),diff(find(phase<0))]))>1
    warning('********* Problem with alternating phase ***********')
end

% $$$ sub_demod(X1,Fc,'Ch4')
% $$$ sub_demod(X2,Fc,'Ch5')
% $$$ sub_demod(X3,Fc,'Ch6')
% $$$ sub_demod(X4,Fc,'Ch7')

% $$$ sub_demod(X1-X2,Fc,'Ch4 - Ch5'), ylim([-.03,.03])


[Tp,mdV,sdV,mEP,mu]=sub_demod(Xinj,Fc,tstr,if_plot, hasTrig, thres, ppp);
         if(if_plot)
                subplot(3,1,2), ylim(y1)
                subplot(3,1,3), ylim(y2)
                saveas(gcf,fname)
         end;
% $$$ sub_demod(X3-X4,Fc,'Ch6 - Ch7'), ylim([-.03,.03])