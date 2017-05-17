function [H,p,ci,stats] = ttest2_printout(sc1,sc2,varargin)
% :Usage:
% ::
%
%     [H,p,ci,stats] = ttest2_printout(sc1,sc2,[doplot],[covts])
%
% one or two sample t-test printout and plot
% covariates are not done yet!
%
% :Inputs:
%
%   **sc1:**
%        data from first group
%
%   **sc2:**
%        data from second group (if missing or empty, performs one-sample
%        t-test)
%
% ..
%    tor wager, last updated Sept 2007 (cosmetic update)
% ..

if length(varargin) > 2
    % covariates of no interest
    X = varargin{2};
    % NOT DONE YET.
end
    
if nargin == 1 || isempty(sc2)
    [H,p,ci,stats] = ttest(sc1,0,.05,'both');
    fprintf(1,'u1 = %3.2f, t(%3.1f) = %3.2f, p = %3.4f\n',nanmean(sc1),stats.df,stats.tstat,p);
else
    [H,p,ci,stats] = ttest2(sc1,sc2,.05,'both','unequal');
    fprintf(1,'u1 = %3.2f, u2 = %3.2f, udiff = %3.2f, t(%3.1f) = %3.2f, p = %3.4f\n',nanmean(sc1),nanmean(sc2),nanmean(sc1) - nanmean(sc2),stats.df,stats.tstat, p);
end

means = [nanmean(sc1) nanmean(sc2)];
stats.means = means;

if length(varargin) > 0 && varargin{1}
    
    hh = bar(means);
    set(hh,'FaceColor',[.7 .7 .7]);
    
    se = stats.sd ./ sqrt(length(sc1)+length(sc2));

    se = se';
    tor_bar_steplot(means,se,{'k'});
    
    set(gca,'XTick',[1 2],'XTickLabel',{'High' 'Low'});
        
end


return

