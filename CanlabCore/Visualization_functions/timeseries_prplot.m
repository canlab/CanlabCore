function [r,bb,f, han] = timeseries_prplot(y,X,cols,varargin)
% :Usage:
% ::
%
% Partial residual plot of a subset [cols] of columns regressing y on X
%
%    [r,bb,f] = timeseries_prplot(y,X,cols,varargin)
%
% - Plots timeseries data (y') removing other regressors (~cols) against partial fitted response (X[cols])
% - Makes sure there is an intercept at the end, whether one was entered or not.
% - This is *not* a partial regression plot
%
% See also: https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/partregr.htm#:~:text=Partial%20regression%20plots%20are%20most,independent%20variables%20in%20the%20model).
% "Partial regression plots are most commonly used to identify leverage
% points and influential data points that might not be leverage points. 
% Partial residual plots are most commonly used to identify the nature of the 
% relationship between Y and Xi (given the effect of the other independent variables 
% in the model)"
%
% :Inputs:
%
%   **y:**
%        y is n x 1 data points, X is n x k model matrix of predictors
%
%        y is adjusted to remove all effects OTHER THAN those columns of X
%        specified in cols, creating a partial residual vector y'.
%
%   **X:**
%        The fitted response to X(:,cols) is plotted against y' to graphically assess the
%        effect of particular columns of X on y.
%
%   varargin includes two optional arguments (enter both *together*):
%   These will create "boxcars" showing trials onsets, e.g., for fMRI
%     1) a vector of trial onsets (to shade in plot)
%     2) length of elements to shade after trial onset.
%
% :Outputs:
%
%   **r:**
%        partial residuals
%
%   **bb:**
%        betas
%
%   **f:**
%        partial fitted response
%
% :Examples:
% ::
%
%    timeseries_prplot(Yfla,X,[2 4],x2,18);
%
% to average across sessions 1 and 2:
%
% ..
%    Tor Wager
% ..

% Make sure there is an intercept at the end, whether one was entered or not.
X = intercept(X, 'end');

% get betas from full model
b = pinv(X) * y; 
bb = b;

% get fit from subset of predictors of interest here

f = X(:,cols) * b(cols);

% get partial residuals, excluding the subset of interest

X(:,cols) = []; b(cols) = [];
r = y - X * b;

% plot
tor_fig; hold on;
han = plot(r,'k','LineWidth', 2); 
han(2) = plot(f,'b','LineWidth', 2);
legend({'Partial residual (data)' 'Partial fit'});

if length(varargin) > 0
    x2 = varargin{1}; len = varargin{2};
    sc = 2 * max(abs(r));
    
    for i = 1:length(x2), hh(i) = fill([x2(i) x2(i)+len x2(i)+len, x2(i)],sc*([0 0 1 1]-.5),[.5 .5 .5]); end
    
    set(hh,'EdgeColor','none','FaceAlpha',.5)
end
    
return

