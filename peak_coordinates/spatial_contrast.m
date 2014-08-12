function [str,dcon] = spatial_contrast(XYZ1,XYZ2)
% function [str,dcon] = spatial_contrast(XYZ1,XYZ2)
%
% this function tests relative locations of individual 
% spatial peaks from clusters.
%
% plots position of XYZ2 relative to XYZ1
% thus, an 'anterior' group position means that XYZ2
% peaks are anterior to XYZ1 peaks
%
% input:  
%	XYZ1,2: n x 3 coordinates (in mm)
%   con:    contrast vector, e.g., [1 -1]
%
% 11/2/03 by Tor Wager
%
% Minor modification/documentation by Tor Wager, Aug 2010
% Now correctly indicates that permutation test is used for p-values.
% See help conf_region for details of the test.


dcon = XYZ2 - XYZ1;
m = mean(dcon); 
mmax = max(abs(max(dcon))) * 2;

mm = mean([XYZ1;XYZ2]);
xx1 = XYZ1 - repmat(mm,size(XYZ1,1),1);
xx2 = XYZ2 - repmat(mm,size(XYZ2,1),1);

figure('Color','w'); subplot(1,2,1);hold on; set(gca,'FontSize',16)
% points
plot3(xx1(:,1),xx1(:,2),xx1(:,3),'rs');
plot3(xx2(:,1),xx2(:,2),xx2(:,3),'go');

hold on; plot3([-mmax mmax],[0 0],[0 0],'k','LineWidth',2); text(-mmax,0,-1,'Left');text(mmax,0,-1,'Right');
hold on; plot3([0 0],[-mmax mmax],[0 0],'k','LineWidth',2); text(0,-mmax,-1,'Posterior');text(0,mmax,-1,'Anterior');
hold on; plot3([0 0],[0 0],[-mmax mmax],'k','LineWidth',2); text(1,0,-mmax,'Inferior');text(1,0,mmax,'Superior');
view(26,34)
xlabel('x (mm)');ylabel('y (mm)'),zlabel('z (mm)')
%camzoom(1.25);
axis vis3d


subplot(1,2,2); hold on; set(gca,'FontSize',16)

% individual ones
for i = 1:size(dcon,1)
    plot3(dcon(i,1),dcon(i,2),dcon(i,3),'b^','MarkerFaceColor','b','MarkerSize',2);
    hh(i) = plot3([0 dcon(i,1)],[0 dcon(i,2)],[0 dcon(i,3)],'b');
end

% confidence sphere
results = confidence_volume(dcon,'r',1);
str=sprintf('%3.0f pts, [%3.1f %3.1f %3.1f], Mean Euclidean Distance from origin = %3.2f, p = %3.4f',size(dcon,1),m(1),m(2),m(3),results.msb, results.pval);
disp(str)

% average vector
hold on; plot3(m(1),m(2),m(3),'b^','MarkerFaceColor','r','MarkerSize',8);
h = plot3([0 m(1)],[0 m(2)],[0 m(3)],'r','LineWidth',4);

% axes
hold on; plot3([-mmax mmax],[0 0],[0 0],'k','LineWidth',2); text(-mmax,0,-1,'Left');text(mmax,0,-1,'Right');
hold on; plot3([0 0],[-mmax mmax],[0 0],'k','LineWidth',2); text(0,-mmax,-1,'Posterior');text(0,mmax,-1,'Anterior');
hold on; plot3([0 0],[0 0],[-mmax mmax],'k','LineWidth',2); text(1,0,-mmax,'Inferior');text(1,0,mmax,'Superior');
view(26,34)
xlabel('x (mm)');ylabel('y (mm)'),zlabel('z (mm)')
%camzoom(1.25);
axis vis3d
title(str)

return
