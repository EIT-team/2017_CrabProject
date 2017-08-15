function [A,phi] = hilbert_demod(x,Fc,Fs)
%
    %tClip = 125*(Fs/2048);
    tClip = 700;
   % Fc=225;
    side = Fc-125;
   % side=50;
   ind = [round((Fc-side)/2),round((Fc+side)/2)];
    
    
    
    
%     if Fc<425
%         ind = [round((Fc-60)/2),round((Fc+60)/2)];
%     else
%         ind = [round(Fc/2-450),round(Fc/2+450)];        
%     end


%      [b,a] = butter(2,[49/Fs,51/Fs],'stop');
%      x = filtfilt(b,a,x);
%      [b,a] = butter(8,(Fc+side)*2/Fs,'low');
%       x = filtfilt(b,a,x);
%      
    
    win = hanning(ind(2)-ind(1)+1);
    
    h = hilbert(x);

    F1 = fft(h);
    F2 = zeros(size(F1));
    F2(ind(1):ind(2),:) = F1(ind(1):ind(2),:).*repmat(win,1,size(x,2));
    h = ifft(F2);

%      [b,a] = butter(4,105*2/Fs,'high');
%      x = filtfilt(b,a,x);

    
    A = abs(h);
%       %%%%%%%%%%
%       [b,a] = butter(6,1200/Fs,'low');
%       A = filtfilt(b,a,A);
%       %%%%%%%%%    
    
    phi = unwrap(angle(h));
% $$$     phi(tClip+1:end-tClip,:) = detrend(phi(tClip+1:end-tClip,:));

    A(1:tClip,:) = repmat(mean(A(tClip+1:end-tClip-1,:),1),tClip,1);
    A(end-tClip+1:end,:) = repmat(A(1,:),tClip,1);

% $$$     phi(1:tClip,:) = repmat(phi(tClip+1,:),tClip,1);
% $$$     phi(end-tClip+1:end,:) = repmat(phi(end-tClip,:),tClip,1);
