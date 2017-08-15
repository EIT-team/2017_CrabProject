function y = plot_EP(Biosemi_fname,ref_chanel, File_path)
%
    Bad_channels = [1:2,13:29];

    % Assume current directory if no File_path
    if nargin==2
        File_path = '.';
    end

    % Fs: sampling frequency
    % T_trig: trigger interval
    Fs = 16384;
    T_trig = Fs/2;
    T_lag = 2458;
    
    Trig_thres = 1e4;
    elec = setdiff(1:29,Bad_channels);
    
    % Load Biosemi data
    EEG = pop_biosig(fullfile(File_path, Biosemi_fname));
    EEG.data=EEG.data-repmat(EEG.data(ref_chanel,:),32,1);
    Data = EEG.data';
    clear EEG
    
    % Data is a NxK matrix
    % where N is the number of samples (recording time)
    % and K is the number of channels 
    % --> select electrode / trigger channels
    Data = Data(:,[1:29,32]);

    % set noisy channels to zero
    Data(:,Bad_channels) = 0;
% $$$     Data = detrend(Data,'constant');

    % This sets a threshold for the trigger channel 
    % NOTE: Might need to adjust threshold, do 'plot(Data(:,end))'
    TrigData = diff(Data(:,end));
    t = find(TrigData>Trig_thres);

    % Gets rid of contiguous time points
    t(find(t(2:end)==(t(1:end-1)+1))+1)=[];

    % Gets rid of the last element (might not have enough recorded
    % data after last trigger)
    t(end)=[];
    t(t<T_lag)=[];

    % Sanity check
    disp(sprintf('min interval=%.2fs, EPs = %i',min(diff(t))/Fs,length(t)))

    % x: electrode channels
    x = Data(:,1:end-1);

    % Remove stimulation artefact
    for i=1:size(x,2)
        for j=1:length(t)
            x(t(j)+(-16:16),i) = mean(x(t(j)+[-17,17],i));
        end
    end

    % Removes linear trend of each channel (DC component)
    x = double(detrend(x,'constant'));

% $$$     % Remove 50-hz mains hum
% $$$     [b1,a1] = butter(3,2*[48,52]/Fs,'stop');
% $$$     x = filtfilt(b1,a1,x);
    
    
% $$$     % Uses least active channel as reference
% $$$     sd_x = std(x(:,elec));
% $$$     [~,i] = min(sd_x);
% $$$     x(:,elec) = x(:,elec)-repmat(x(:,14),1,length(elec));
    
    % Find instances in which Trigger failed
    trig_Fail = nnz(diff(t)>1.1*T_trig);
    if any(trig_Fail)
        disp(sprintf('Trigger failed %i times!',trig_Fail))
    end
    
    % y: average EP (from 0 to 0.25s after trigger)
    nElec = size(x,2);
    nRep = length(t);
    X = cell(nElec,1);
    y = zeros(T_trig,nElec);
    for i=setdiff(1:nElec,Bad_channels)
        x_i = zeros(T_trig,nRep);
        for j=1:nRep
            x_i(:,j) = x(t(j)-T_lag+(1:T_trig),i);
        end
        x_i = x_i-repmat(x_i(T_lag+1,:),T_trig,1);
        y(:,i) = mean(x_i,2);
        X{i} = x_i;
    end
    y(1:T_lag,:) = [];


    maxV = max(reshape(y(1:800,:),[],1));
    [I,J] = ind2sub([800,size(y,2)],find(y(1:800,:)==maxV));

    disp(sprintf('%s: Max. EP = %.3f mV (electrode #%i, time=%.1f ms)',...
                 Biosemi_fname,maxV/1e3,J,I/Fs*1e3));

    figure('Position',get(0,'ScreenSize'));

    % plot the first ~50ms of all channels (y: average EP)
    subplot(2,1,1);
    plot((1:800)/Fs*1e3,y(1:800,:));
    xlabel('Time (ms)');
    ylabel('Average EP (uV)');
    title(sprintf('%s: Max. EP = %.3f mV (electrode #%i, time=%.1f ms)',...
                  Biosemi_fname,maxV/1e3,J,I/Fs*1e3),'Interpreter','none');

    % image representation of average EPs
    subplot(2,1,2);
    imagesc(y(1:800,:)');
    colorbar;
    xlabel('Time (bins)');
    ylabel('Electrode #');

    %saveas(gcf,sprintf('EP_%s.png',Biosemi_fname))
