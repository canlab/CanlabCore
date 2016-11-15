function vif = getvif(X, no_add_intercept, varargin)

% function vif = getvif(model design matrix (X), [no_add_intercept], varargin)
%
% Optional arguments:
%
% 1.  first varargin must be flag for no_add_intercept.
% The only case where we may not want to add an intercept is if we already
% have intercepts for each subset of observations (e.g., each run in fmri)
% and these are sufficient.
%
% 2.
% pass in 'wh' followed by vector of indices to return the vifs for only those columns,
% though vif will still be calculated on entire model.
%
% refactored and updated documentation: by tor wager, aug 2015
%
% Examples:
% Generate X:
% [X, delta, delta_hires, hrf] = onsets2fmridesign(ons, TR, scanLength, 'hrf', 'parametric_standard', pm);
% OR
% X = mvnrnd([1 1 1 1], eye(4), 100);
%
% vifs = getvif(X)
%
% getvif(X, 0, 'plot')
%
%
% See also:
% scn_spm_design_check

% Remove intercept, we will add later to each model
X = intercept(X, 'remove');         % remove if present

[n, k] = size(X);

% Optional arguments
% -------------------------------------------------------------------------

if nargin < 2, no_add_intercept = 0; end

ind = strcmp('wh', varargin);

if ~isempty(ind) && sum(ind) ~= 0
    
    wh=varargin{ind+1};
    if size(wh, 1) > size(wh, 2), wh=wh';end %transpose if needed
    
else
    
    wh = 1:k;
    
end

doplot = 0;

if any(strcmp('plot', varargin))
    doplot = 1;
end

% calculate
% -------------------------------------------------------------------------

% initialize, exclude the intercept
vif = [];  %ones(1, length(wh));

intcpt = ones(n, 1);

for i = wh
    
    if no_add_intercept
        Xi = X;
    else
        Xi = [X intcpt];
    end
    
    % check that model has an intercept, or that all preds are
    % mean-centered.  if both of these conditions are false, Rsquared is uninterpretable see [http://statistiksoftware.blogspot.com/2013/01/why-we-need-intercept.html]
    % and thus VIF is uninterpretable.  Yoni, May 2015
    hasint = any(var(Xi)==0);
    totallymc = sum(mean(Xi))==0;
    
    if ~hasint && ~totallymc
        warning('Model has no intercept, and model is not mean-centered; VIF is not interpertable');
        vif = NaN * zeros(1, length(wh));
        return
    end
    
    
    y = Xi(:, i);
    
    Xi(:, i) = [];
    
    b = Xi\y;
    fits = Xi * b;
    
    rsquare = var(fits) / var(y);
    
    % sometimes rounding error makes rsquare>1
    if rsquare >= 1, rsquare = .9999999; end
    
    vif(end + 1) = 1 / (1 - rsquare);
    
    
end  % regressor

if doplot
    
    hold on;
    plot(wh, vif, 'ko', 'MarkerFaceColor', [1 .5 0], 'MarkerSize', 8);
    
    ylabel('Variance inflation factor'); xlabel('Predictor number');
    plot_horizontal_line(1, 'k');
    h = plot_horizontal_line(2, '--'); set(h, 'Color',[0 .3 .7]);
    h = plot_horizontal_line(4, '--'); set(h, 'Color', [0 .7 .3]);
    h = plot_horizontal_line(8, '-'); set(h, 'Color', [.7 .3 0]);
    
    disp('Variance inflation: 1 (black line) = minimum possible (best)');
    disp('Successive lines indicate doublings of variance inflation factor.');
    disp('Red boxes have extremely high VIFs, perfect multicolinearity');
    
    title('Var. Inflation (VIFs) in full design');
    
    mymax = max(10, max(vif(vif < 1000))+3);
    
    set(gca, 'YLim', [0 mymax], 'FontSize', 20, 'XTick', 1:length(wh), 'XLim', [min(wh)-.5 max(wh) + .5]);
    
    whhigh = find(vif > 1000);
    for i = 1:length(whhigh)
        
        drawbox(wh(whhigh(i)) - .3, .6, 0, mymax, [1 0 0]);
        
    end
    
end

end % function
