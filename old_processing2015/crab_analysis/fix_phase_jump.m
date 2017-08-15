function Data = fix_phase_jump(Data,info)
%

    TrigData = detrend(Data(:,info.Trig_Chan),'constant');
    TrigData(1:13) = 0;
    TrigData(end-13:end) = 0;

    Ns = floor(size(TrigData,1)/info.Fs);% number of integer seconds
    Time = [1:length(TrigData)]/info.Fs;
    
    info.Th_Trig = 1e3;%median(max((reshape(TrigData(1:Ns*info.Fs),info.Fs,Ns))))/2;
    info.Frame = ceil(info.Fs/(info.Stim_Rate)); % for SD Card

    ff = find(TrigData>info.Th_Trig);
    
    % the following line was added to reduce 5 close pulses
    % appearing in the AER trig into 1 pulse
    ff(find(diff(ff)<30)+1) = [];

    % the following blocks corrects for the differnt analog properties
    % of the trigger signal coming from teh current source.
    % it is stonger and wider and therefore the peak might be in the
    % folloing 1-2 samples.    
    for k=1:length(ff)
        [tmp1,ind] = max(TrigData(ff(k)-12:ff(k)+12));
        ff(k) = ff(k)+ind-13;
    end
    
    Phase_tmp = rem(ff,info.Frame);
    info.Phase = Phase_tmp(1);
    pj_ind_plus = find(diff(Phase_tmp)==1); % extra sample
    pj_ind_minus = find(diff(Phase_tmp)==-1); % missing sample

    if ~isempty(pj_ind_plus),
        disp('+')
        TrigData(ff(pj_ind_plus+1)) = [];
        Data(ff(pj_ind_plus+1),:) = [];    
    end
    if ~isempty(pj_ind_minus),
        disp('-')
        Phase_Jump_Location_Minus = ff(pj_ind_minus)-2457;  %%%%%%%%%%%%%%%%%%%
        for r = 1:length(Phase_Jump_Location_Minus)
            TrigData = [TrigData(1:Phase_Jump_Location_Minus(r)+(r-1)); ...
                TrigData(Phase_Jump_Location_Minus(r)+1+(r-1)); ...
                TrigData(Phase_Jump_Location_Minus(r)+1+(r-1):end)];
            Data = [Data(1:Phase_Jump_Location_Minus(r)+(r-1),:); ...
                Data(Phase_Jump_Location_Minus(r)+1+(r-1),:); ...
                Data(Phase_Jump_Location_Minus(r)+1+(r-1):end,:)];
        end
    end

    %%%
    Trig_Delay=150;

    Phase_A = floor(info.Phase-Trig_Delay/1000*info.Fs);
    if Phase_A<0, Phase_A=Phase_A+info.Frame;end
    Ind_A = 1+Phase_A:length(Data);
    N_A = floor(length(Ind_A)/info.Frame);
    if N_A/2~=floor(N_A/2), N_A=N_A-1;end

    Data = detrend(double(Data(Ind_A(1:N_A*info.Frame),:)));
    