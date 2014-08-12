function y2 = discrim_plot(discf,y,varargin)
% y2 = discrim_plot(discf,y,[new figure],[column to do partial corr plot of])
%
% Plots correlation or partial correlation
% Color codes by median split of y
%
% see also:
% cluster_discrim
% cluster_discrim_montage

dofig = 1; dopr = 0;
if length(varargin) > 0, dofig = varargin{1};,end
if length(varargin) > 1, dopr = varargin{2};,end

y2 = mediansplit(y);

if dopr
    [discf(:,1),y,r,p,rrob,prob] = partialcor(discf,y,dopr);
    fprintf(1,'\nCalculating partial correlation.\n')
end


if dofig, figure('Color','w');,end

plot_correlation_samefig(discf(:,1),y,[],'ko',0,1);
wh = find(y2>0); wh2 = find(y2<0);
hold on; plot(discf(wh2,1),y(wh2),'bo','MarkerSize',10,'LineWidth',2);
hold on; plot(discf(wh,1),y(wh),'ro','MarkerSize',10,'LineWidth',2);
xlabel('Discriminant function'); ylabel('Behavior');

return
