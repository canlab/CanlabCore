function [x,y,r,p,se,meany,stats] = partialcor(X,y,i,doprint,dorobust)
% function [x,y,r,p,se,meany,stats] = partialcor(X,y,i,[doprint],[dorobust])
%
% partial correlation between column i of X and y
% X is n x k matrix of predictors (you can include an intercept; will be added if missing)
% y is n x 1 data vector
%
% x, y are adjusted predictor and data
% r, p are Pearson's correlation and p-value
% ***** Default option: use robust IRLS. ******
%
% Warning: p-values are not adjusted for use of degrees of freedom in
% partial correlation computation w/o robust option!  OK for robust
%
% Modified by Tor Wager, July 2006
% :mean of y retained in case we're interested in the intercept
% :if X of interest is intercept, returns robust mean, controlling for
% other columns of X
% if robust outputs are requested, y is created from robust stats

if nargin < 5, dorobust = 1; end
if nargin < 4, doprint = 0; end

fullX = X;
yorig = y;

% ---------------------------------------------------
% Prepare and check X matrix
% ---------------------------------------------------

% x and y are adjusted predictor and data
x = X(:,i);
xorig = x;

% add intercept if not already in model
no_intercept = ~any(all(diff(fullX) < eps));
if no_intercept, fullX(:,end+1) = 1; end

wh_intercept = find(all(abs(diff(fullX)) < eps)); %8/5/15 added abs() - Scott Schafer
%abs() is necessary to prevent mis-identifying between group contrasts as
%an intercept (because going from 1 -> -1 makes diff = -2, which is "less
%than" zero.

% center predictors (except intercept)
% not needed; and this changes the interpretation of the betas in case x is
% not centered.
%fullX = scale(fullX,1);
%fullX(:,wh_intercept) = 1;

redX = fullX;
redX(:,i) = [];     % reduced model without param(s) of interest

if ~(wh_intercept == i)
    % we are removing the intercept, so we want to add it in later
    wh_intercept_red = find(all(abs(diff(redX)) < eps)); %8/5/15 added abs() - Scott Schafer
end


% ---------------------------------------------------
% Fit GLM model to get stats
% ---------------------------------------------------
if ~dorobust
    % OLS
    rstr = 'OLS ';
    [b,dev,stats] = glmfit(fullX,y,'normal','constant','off');
    stats.w = ones(size(y,1),1);
    stats.b = b;

else
    % robust IRLS
    rstr = 'Robust ';

    [b,stats]=robustfit(fullX,y,'bisquare',[],'off');
    stats.b = b;
end

% ---------------------------------------------------
% % Adjust x and y for regressors of no interest
% ---------------------------------------------------

% weighted fit-forming matrix for no-interest cols of X
[bb,invxwx,bform,fits] = get_betas_singleweight(redX,y,stats.w);
y = y-fits;

bbx =  bform * x;
x = x - redX * bbx;

if ~(wh_intercept == i)
    % we are removing the intercept, so we want to add it in to show mean
    y = y + bb(wh_intercept_red);
    x = x + bbx(wh_intercept_red);
end

% this code gives the same results,but not for adding in intercept
% W = diag(stats.w);
% %%%% this is slightly different: fitmtx = W * redX * pinv((W * redX));
% fitmtx = redX * pinv((W * redX));
% x = x - fitmtx * x + %%%but not here: must add rob mean mean(x);       % keep mean b/c it's of interest in plot
% y = y - fitmtx * y + mean(y); 

if wh_intercept == i
    % intercept value
    r = stats.b(i);
else
    r = partial_r_wt(x,y,stats.w);
end

% ---------------------------------------------------
% % Opt: Model comparison version
% something seems wrong with this.
% ---------------------------------------------------
% Model comparison method
domodelc = 0;
if domodelc
    if ~isempty(redX)
        [r,r2,stats] = partial_r(fullX,redX,stats,yorig,wh_intercept,i,dorobust);
        y = stats.resid + x * b(i);     % adjusted y
    else
        % undefined, we have only the intercept
        r = 0;
    end
end


% this is the mean estimate because predictors are centered
meany = stats.b(wh_intercept);
t = stats.t(i);
p = stats.p(i);
se = stats.se(i);

% very low P-values: do not use exactly zero, because of problems
p = max(p, 10 * eps);

if wh_intercept == i  %all(x - mean(x) < eps)
    % x is an intercept; weighted mean
    %r = stats.w' * y ./ sum(stats.w);
    pstr = 'intercept val';
else
    % weighted correlation
    %r = weighted_corrcoef([x y],stats.w);
    %r = r(1,2);
    pstr = 'partial r';
end

if doprint
    fprintf(1,'%s %s: %3.2f, b = %3.2f, se = %3.2f, t(%3.0f) = %3.2f, p = %3.4f, d = %3.2f\n',rstr,pstr,r,b(i),se,stats.dfe,stats.t(i),p, stats.t(i) ./ (size(X, 1) .^ .5));
end



return



function r = partial_r_wt(x,y,w)

r = weighted_corrcoef([x y],w);
r = r(1,2);

return






function [r_reduced,r2_reduced,stats] = partial_r(fullX,redX,stats,y,wh_intercept,wh_interest,dorobust)

% "Model comparison" version
% pass in yorig

if ~dorobust
    % OLS
    rstr = 'OLS ';
%     [b,dev,stats] = glmfit(fullX,y,'normal','constant','off');
%     stats.w = ones(size(y,1),1);
%     stats.b = b;
    
    if ~isempty(redX)
        [br,dev,statsr] = glmfit(redX,y,'normal','constant','off');
    end
    statsr.w = ones(size(y,1),1);

else
    % robust IRLS
    rstr = 'Robust ';

%     [b,stats]=robustfit(fullX,y,'bisquare',[],'off');
%     stats.b = b;
    
    if ~isempty(redX)
        [br,statsr]=robustfit(redX,y,'bisquare',[],'off');
    end
end

mysign = sign(stats.b(wh_interest));

y = y - stats.b(wh_intercept);  %mean(y);    % deviations; only OK of covariates are centered!

W = diag(stats.w);
W = W' * W;         % squared weights

SSt = y' * W * y;   % total variability in data

e = stats.resid;
er = statsr.resid;

SSf = e' * W * e;   % full model sum sq. error

SSr = er' * W * er;  % reduced model sse

r2_reduced = (SSr - SSf) ./ SSt;    % R2 accounted for by predictor of interest

r_reduced = sqrt(r2_reduced);       % partial correlation (absolute value)

r_reduced = r_reduced * mysign;

%r2_full = (SSt - SSf) ./ SSt;       % full model r-square

%r_full = sqrt(r2_full);              % multiple R for full regression

return







function [betas,invxwx,bform,fits] = get_betas_singleweight(X,Y,w)

W = diag(w);                    % Weight matrix

%X = repmat(1,m,1);              % Design matrix - 1 column of all ones to calculate average
% and, separately, use bcon if that's entered

invxwx = inv(X'*W*X);
bform = invxwx * X'* W;         % beta-forming matrix.  hat = X * bform

% rows are columns of design (X), cols are Y variables
betas = bform*Y;

if nargout > 3
    fits = X * betas;
end

return

% if ~dorobust
%     % X is no interest; x is of interest
%     X(:,i) = [];
%
%     W = X * pinv(X);
%
%     x = x - W * x;
%
%     % leave mean of y in in case we're interested in the intercept
%     y = y - W * y + mean(y);


