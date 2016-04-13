function cl = iimg_princomp_display(volInfo,k,overlay,dofigs)
% :Usage:
% ::
%
%     cl = iimg_princomp_display(volInfo,k,overlay,dofigs)
%
% Show output from iimg_princomp
%
% :Example:
% ::
%
%    cl2 = iimg_princomp_display(volInfo2,1,EXPT.overlay);



tor_fig; plot(volInfo.eigval(:,k),'ko-','LineWidth',2);

%overlay = [];
%k = 1;

title('Eigenvalues');

indx = double(volInfo.image_indx);
mycorrs = volInfo.corr_with_comps(:,k);
mycorrs(abs(mycorrs) < .5) = 0;
indx(volInfo.wh_inmask) = mycorrs;

voldata = iimg_reconstruct_3dvol(indx,volInfo);
cl = mask2clusters(double(voldata),volInfo.mat);

if dofigs
montage_clusters(overlay,cl,[2 2]);
%montage_clusters(overlay,cl,{'y'});

surfh = cluster_surf(cl,5,'heatmap');
lightRestoreSingle(gca); axis off; set(gcf,'Color','w'); material dull
end


%% in progress: make scatterplot with behaviora
% tor_fig; r = prplot(beh,scale(X),1,1,{'yo'});
% xlabel('Average activity in network (z-score)');
% ylabel('Reappraisal success (reported)');

