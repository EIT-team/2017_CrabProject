electrodes = [5 9 13 17 5];
figure
for i = 1:5
    dist_left = (electrodes(i)-2)*4;
    dist_right = (electrodes(i)-1)*4;
    av_dist = (dist_left+dist_right)/2
    dVm = mean(dVstruc(i).dVp,2);
    subplot(5,1,i)
    hold on
    plot(T,dVstruc(i).dVp,'color',[0.7 0.7 0.7]);
    plot(T,dVm,'linewidth',3);
    title(sprintf('dV at: %d mm\n',av_dist));
    ylabel('uV');
    xlabel('T ms');
    hold off
    xlim([-5 100]);
    ylim([-0.2 0.2]);
end