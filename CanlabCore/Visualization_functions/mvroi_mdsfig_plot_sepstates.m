function mvroi_mdsfig_plot_sepstates(CLUSTER,CORRELS,SPEC)
% ::
%
%    mvroi_mdsfig_plot_sepstates(CLUSTER,CORRELS,SPEC,DATA.CORRELS.AVGSTATS.sigmat_uncorrected,DATA.CORRELS.DIFSTATS.sigmat_uncorrected,titlestr)
%
% Plots SEPARATE STATES
%
% :Example: 
% ::
%
%    mvroi_mdsfig_plot2(DATA.CLUSTER,DATA.SPEC,DATA.CORRELS.AVGSTATS.sigmat_uncorrected,DATA.CORRELS.DIFSTATS.sigmat_uncorrected,'Uncorrected');
%
scnsize = get(0,'ScreenSize');
figure('position',[50 50 scnsize(3)-100 scnsize(4)/2],'color','white')

nstates = SPEC.numstates;

for i = 1:nstates
    
    subplot(1,nstates,i);
    set(gca,'FontSize',14)
    
    titlestr = SPEC.tasknames{i};
    
    su = CORRELS.STATESTATS(i).sigmat_uncorrected;
    sc = CORRELS.STATESTATS(i).sigmat;
    
    mdsfig(CLUSTER.Gs(:,1:2),CLUSTER.names,CLUSTER.classes,su); 
    title(['State ' num2str(i) ', ' titlestr ' (p < .05 uncorrected)'])
end
drawnow

figure('position',[50 50 scnsize(3)-100 scnsize(4)/2],'color','white')

for i = 1:nstates
    
    subplot(1,nstates,i);
    set(gca,'FontSize',14)
    
    titlestr = SPEC.tasknames{i};
    
    su = CORRELS.STATESTATS(i).sigmat_uncorrected;
    sc = CORRELS.STATESTATS(i).sigmat;
    
    mdsfig(CLUSTER.Gs(:,1:2),CLUSTER.names,CLUSTER.classes,sc); 
    title(['State ' num2str(i) ', ' titlestr ' (p < .05 corrected)'])
end

drawnow

return

