%%
N1=load('E:\crabcheck\2017-11-21 - Crab Again\EIT_001_out.mat'); %worked
N2=load('E:\crabcheck\2017-11-21 - Crab Again\EIT_squirt_001_out.mat'); %worked again
N3=load('E:\crabcheck\2017-11-21 - Crab Again\171122-1-N5_EIT_002_out.mat'); %didnt
N4=load('E:\crabcheck\2017-11-21 - Crab Again\171122-12-N6_EIT_003_1ma_out.mat'); %no change at all
N5=load('E:\crabcheck\2017-11-21 - Crab Again\bigc_out.mat'); %positive and short nerve

labels={'Working','Working2','NoChange New','No change','20ua short positive'};
xlims = [-20 40];

good_elec=1;


%% Boundary voltages

figure
hold on
plot(N1.EIT_BVm);
plot(N2.EIT_BVm);
plot(N3.EIT_BVm);
plot(N4.EIT_BVm);
plot(N5.EIT_BVm);
hold off 
xlabel('Electrode');
ylabel('BV uV')
title('Boundary Voltages')
legend(labels);

%%
figure
hold on
plot(N1.Out(good_elec).T,N1.Out(good_elec).EPm);
plot(N2.Out(good_elec).T,N2.Out(good_elec).EPm);
plot(N3.Out(good_elec).T,N3.Out(good_elec).EPm);
plot(N4.Out(good_elec).T,N4.Out(good_elec).EPm);
plot(N5.Out(good_elec).T,N5.Out(good_elec).EPm);
hold off
xlim(xlims)
title('EP')
ylabel('uV')
xlabel('T ms');
legend(labels);

%%

figure
hold on
plot(N1.Out(good_elec).T,N1.Out(good_elec).dVm);
plot(N2.Out(good_elec).T,N2.Out(good_elec).dVm);
plot(N3.Out(good_elec).T,N3.Out(good_elec).dVm);
plot(N4.Out(good_elec).T,N4.Out(good_elec).dVm);
plot(N5.Out(good_elec).T,N5.Out(good_elec).dVm);
hold off
xlim(xlims)
title('dV')
ylabel('uV')
xlabel('T ms');
legend(labels);


%%

figure
hold on
plot(N1.Out(good_elec).T,mean(N1.Out(good_elec).dV_sig_orig,2));
plot(N1.Out(good_elec).T,mean(N2.Out(good_elec).dV_sig_orig,2));
plot(N1.Out(good_elec).T,mean(N3.Out(good_elec).dV_sig_orig,2));
plot(N1.Out(good_elec).T,mean(N4.Out(good_elec).dV_sig_orig,2));
plot(N1.Out(good_elec).T,mean(N5.Out(good_elec).dV_sig_orig,2));
hold off
xlim(xlims)
title('dV')
ylabel('uV')
xlabel('T ms');
legend(labels);