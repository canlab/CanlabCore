function [r,bb,f] = timeseries_prplot(y,X,cols,varargin)
% :Usage:
% ::
%
%    [r,bb,f] = timeseries_prplot(y,X,cols,varargin)
%
% Plots timeseries data (y') against fitted response (X)
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
%   varargin includes two optional arguments:
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



b = pinv(X) * y; % get betas
bb = b;

% get fits

f = X(:,cols) * b(cols);

% get partial residuals

X(:,cols) = []; b(cols) = [];
r = y - X * b;

% plot
tor_fig;
plot(r,'k','LineWidth',2); hold on;
plot(f,'k--');

if length(varargin) > 0
    x2 = varargin{1}; len = varargin{2};
    sc = 2 * max(abs(r));
    for i = 1:length(x2), hh(i) = fill([x2(i) x2(i)+len x2(i)+len, x2(i)],sc*([0 0 1 1]-.5),[.5 .5 .5]);,end
    set(hh,'EdgeColor','none','FaceAlpha',.5)
end
    
return

