function nmdsfig_legend(X,r)
% X is stim coords, r is correlation coefficient matrix
% nmds figure should be current fig.
%
% nmdsfig_legend(c.ClusterSolution.X,c.r)

% get x limit and dims of current nmds fig.
xsz = get(gcf,'Position'); xfig = xsz(1); xsz = xsz(3); 
xlim = get(gca,'XLim');
xlim = [-.05 diff(xlim)-.05];  % shift it over

d=pdist(X)'; 
r=r-eye(size(r)); r=squareform(r)';
b = pinv([r ones(size(r,1),1)]) * d;

x = [.9 .7 .5]; %0 -.5 -.7 -.9];
n = length(x);
y = b(1) * x + b(2);

figure('Color','w');
hold on;

shiftx = .05;

for i = 1:n
    if x(i) < 0, mycol = 'c'; else mycol = 'k'; end
    
    plot([0 y(i)],[i i],mycol,'LineWidth',2);
    text(y(i) + shiftx,i,['r = ' num2str(x(i))],'FontSize',14)
end

po = get(gcf,'Position');
po(1) = xfig;
po(3) = xsz;
po(4) = 80; % fix yaxis to be small
set(gcf,'Position',po);
set(gca,'XLim',xlim,'YLim',[-.5 n+.5]);

axis off
scn_export_papersetup(100)

return
