function STEPWISE = stepwise_tor(dat,y,varargin)
% :Usage:
% ::
%
%     STEPWISE = stepwise_tor(dat, y, [pred. names], [alpha])
% 
% Stepwise regression using Matlab with a couple of extras:
%
% Print omnibus F-values for stepwise regression
%
% get adjusted R-squared
%
% save output structure
%
% ..
%    tor wager
% ..

alph = .05;
if length(varargin) > 1
    alph = varargin{2};
end
    
if length(varargin) > 0
    nms = varargin{1};
else
    for i = 1:size(dat,2), nms{i} = ['V ' num2str(i)]; end
end

    % Add logistic regression for categorical 1/0 DVs!
    

    [STEPWISE.b,STEPWISE.se,STEPWISE.pval, ...
    STEPWISE.inmodel,stats] = ...
    stepwisefit(dat,y,'penter',alph,'display','off');

    t = stats.TSTAT;
    STEPWISE.t = t;
    
    S = STEPWISE;
    tabled = [S.b S.se t S.inmodel' S.pval]; tabnames = {'Beta' 'Std. Err.' 't-value' 'In Model' 'p-value'};
    print_matrix(tabled, tabnames,nms);
    
    N = stats.dfe+stats.df0;
    r2 = (stats.SStotal-stats.SSresid) ./ stats.SStotal;
    adjr2 = 1 - (1-r2)*((N - 1) ./ (N - stats.df0 - 1));
    
    fprintf(1,'\nOmnibus F(%3.0f,%3.0f) = %3.2f, RMSE = %3.2f, p = %3.6f, Adj R^2 = %3.2f\n', ...
    stats.df0,stats.dfe,stats.fstat, ...
    stats.rmse,stats.pval,adjr2);
    stats.adjr2 = adjr2;
    STEPWISE.stats = stats;
    
    

    
 return
    
