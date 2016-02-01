function stats = loess_multilevel(X, Y, varargin)
% :Usage:
% ::
%
%     stats = loess_multilevel(X, Y, varargin)
%
% Options:
%  - case 'fit_and_average', meth = 'fit_and_average';
%  - case 'regularization', regularization = varargin{i+1};
%  - case 'order', loessorder = varargin{i+1};
%  - case {'robust', 'dorobust'}, dorobust = varargin{i+1};
%  - case {'nboot', 'bootsamples'}, nboot = varargin{i+1};
%  - case 'plot', doplot = 1;
%  - case {'plotall'}, plotsummary = 0; doplot = 1;
%  - case 'color', color = varargin{i+1};
%  - case {'existingfig', 'samefig'}, newfig = 0;
%
% :Example:
% ::
%
%    stats = loess_multilevel(models{i}.X, models{i}.Y, 'regularization', .8, 'order', 1, 'color', [.8 .3 0], 'plot', 'samefig');
%
% :See Also: scn_stats_helper_functions

if ~exist('loess.m', 'file')
    disp('You must have loess.m from the Dataviz toolbox on your path.');
    return
end

N = length(X);

regularization = .8;
loessorder = 1; % 1 or 2, linear or quadratic
dorobust = 0; % 0 or 1, bisquare weighting
doplot = 0;
plotsummary = 1;
color = 'b';
newfig = 1;
meth = 'bootstrap'; % pool and bootstrap; alternative: fit individually and average, 'fit_and_average'
nboot = 30;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            % functional commands
            case 'fit_and_average', meth = 'fit_and_average';
            case 'regularization', regularization = varargin{i+1};
            case 'order', loessorder = varargin{i+1};
            case {'robust', 'dorobust'}, dorobust = varargin{i+1};
            case {'nboot', 'bootsamples'}, nboot = varargin{i+1};
            case 'plot', doplot = 1;
            case {'plotall'}, plotsummary = 0; doplot = 1;
            case 'color', color = varargin{i+1};
            case {'existingfig', 'samefig'}, newfig = 0;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


if doplot
    if newfig, create_figure('X - Y plot'); end
    
    if ~plotsummary
        for i = 1:N
            hh(i) = plot(X{i}, Y{i}, 'ko', 'MarkerSize', 4);
        end
    end
    
end

% Set range
% -------------------------------------
allx = cat(1, X{:});
allx = allx(:);
x = [min(allx) max(allx)];
clear allx
xx = linspace(x(1), x(2), 100);

% Get Loess fits
% -------------------------------------
switch meth
    case 'bootstrap'
        X = double(cat(1, X{:}));
        Y = double(cat(1, Y{:}));
        
        bootsam = setup_boot_samples(X, nboot);
        
        %matlabpool open;  % only works with the Parallel Computing Toolbox
        opt = statset('UseParallel','always');
        
        fhan = @(x, y) loess(x, y,  xx, regularization, loessorder, dorobust);
        
        fprintf('Bootstrapping %3.0f samples. ', nboot);
        t1 = tic;
        yy = bootstrp_havesamples(bootsam, fhan, X, Y);
        toc(t1);
        
    case 'fit_and_average'
        for i = 1:N
            
            yy(i, :) = loess(X{i}, Y{i}, xx, regularization, loessorder, dorobust);
            
            % avoid crazy extrapolation -- limit prediction to range of data
            wh_out = xx < min(X{i}) | xx > max(X{i});
            yy(i, wh_out) = NaN;
            
            if doplot && ~plotsummary
                plot(xx, yy(i,:))
            end
            
        end
end


stats.fit_mean = nanmean(yy)';

switch meth
    case 'bootstrap'
        stats.fit_se = std(yy)'; % std of bootstrap dist is ste
    case 'fit_and_average'
        stats.fit_se = ste(yy)';
end

stats.fit_sd = nanstd(yy)';
stats.fit_indiv = yy;
stats.xvals = xx;

if doplot
    stats.error_fill_handle = fill_around_line(stats.fit_mean, stats.fit_se, color, stats.xvals);
    stats.line_handle = plot(stats.xvals, stats.fit_mean, 'Color', color, 'LineWidth', 3);
    stats.error_fill_handle = plot(stats.xvals, stats.fit_mean+stats.fit_se, 'Color', color, 'LineWidth', .5);
    stats.error_fill_handle = [stats.error_fill_handle plot(stats.xvals, stats.fit_mean-stats.fit_se, 'Color', color, 'LineWidth', .5)];
end

end
