function mvroi_mdsfig_plugin2(CLUSTER,SPEC,sigmatavg,sigmatdif)
% ::
%
%    mvroi_mdsfig_plugint2(CLUSTER,SPEC,DATA.CORRELS.AVGSTATS.sigmat_uncorrected,DATA.CORRELS.DIFSTATS.sigmat_uncorrected)
%

scnsize = get(0,'ScreenSize');
figure('position',[50 50 scnsize(3)-100 scnsize(4)/2],'color','white')

subplot(1,2,1);
set(gca,'FontSize',14)
mdsfig(DATA.CLUSTER.Gs(:,1:2),DATA.CLUSTER.names,DATA.CLUSTER.classes,sigmatavg); 
title('Average across conditions')

subplot(1,2,2);
set(gca,'FontSize',14)
mdsfig(DATA.CLUSTER.Gs(:,1:2),DATA.CLUSTER.names,DATA.CLUSTER.classes,sigmatdif); 
title(sprintf(['Contrast [' repmat('%3.0f ',1,size(DATA.SPEC.comps,2)) '] across conditions, %s'],DATA.SPEC.comps(1,:),DATA.SPEC.comptitle{1}))
