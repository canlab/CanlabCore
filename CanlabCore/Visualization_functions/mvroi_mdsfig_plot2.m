function mvroi_mdsfig_plot2(CLUSTER,SPEC,sigmatavg,sigmatdif,titlestr)
% ::
%
%    mvroi_mdsfig_plot2(CLUSTER,SPEC,DATA.CORRELS.AVGSTATS.sigmat_uncorrected,DATA.CORRELS.DIFSTATS.sigmat_uncorrected,titlestr)
%
% :Example:
% ::
%
%    mvroi_mdsfig_plot2(DATA.CLUSTER,DATA.SPEC,DATA.CORRELS.AVGSTATS.sigmat_uncorrected,DATA.CORRELS.DIFSTATS.sigmat_uncorrected,'Uncorrected');
%
scnsize = get(0,'ScreenSize');
figure('position',[50 50 scnsize(3)-100 scnsize(4)/2],'color','white')

subplot(1,2,1);
set(gca,'FontSize',14)
mdsfig(CLUSTER.Gs(:,1:2),CLUSTER.names,CLUSTER.classes,sigmatavg); 
title(['Average across conditions, ', titlestr])

subplot(1,2,2);
set(gca,'FontSize',14)
mdsfig(CLUSTER.Gs(:,1:2),CLUSTER.names,CLUSTER.classes,sigmatdif); 
title(sprintf(['Contrast [' repmat('%3.0f ',1,size(SPEC.comps,2)) '] across conditions, %s'],SPEC.comps(1,:),SPEC.comptitle))


return

