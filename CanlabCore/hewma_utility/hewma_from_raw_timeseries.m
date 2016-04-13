function STATS = hewma_from_raw_timeseries(raw_data, varargin)
% :Usage:
% ::
%
%     STATS = hewma_from_raw_timeseries(raw_data, varargin)
%
% :Inputs:
%
%   **obj:**
%        imraw_data, a N subjects x T time points matrix of data
%
% :Optional Inputs:
%
%   For all optional inputs, enter a keyword followed by an input value/argument
%   (this is a standard Matlab format for entering optional arguments)
%   e.g., 'lam', .3 to enter a lambda value of .3
%
%   Here are the optional inputs and their defaults:
%
%   **lam = .2;**
%        ewma smoothing param
%
%   **L = 2;**
%        ewma control limit
%
%   **noisemodel = 'AR(2)';**
%        noise structure type
%
%   **base_timepts = 60;**
%        baseline time points to use
%
%   **doplot = 1;**
%        plot toggle (1/0)
%
%   **dodetrend = 1;**
%        hewma linear detrending toggle (1/0)
%
%   **grpcontrast = [];**
%        vector of 1 or -1 values for group assignment for
%        each subject
%
%   **samprate = .5;**
%        sampling rate in Hz, 1/TR;  TR = 2 by default
%
% :Note: see ewma5.m and hewma2.m for more information on these inputs
%
% ..
%    Tor Wager, Dec 2008
% ..

lam = .2;
% -----------------------------------------------------------------------
% defaults
% -----------------------------------------------------------------------

L = 2;
noisemodel = 'AR(2)';
base_timepts = 60;
doplot = 1;
dodetrend = 1;
grpcontrast = [];
samprate = .5;  % TR = 2

% -----------------------------------------------------------------------
% optional inputs
% -----------------------------------------------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        % name variables with each input name
        eval([varargin{i} ' = ' varargin{i + 1}]);
        varargin{i + 1} = [];
        
        switch varargin{i}
            case 'rows', rowsz = varargin{i+1};
            %otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% -----------------------------------------------------------------------
% ewma for each subject
% -----------------------------------------------------------------------
[sigmap,sigthresh,chptmap,cntmap,widmap,t1pts,t2pts,t3pts,normdat,stats_ewma] = ...
    ewma5(raw_data, lam, L, noisemodel, base_timepts, doplot);

if doplot
    create_figure('Timeseries', 3, 1); plot_error(raw_data);
    axis tight; title('Raw data');
    subplot(3, 1, 2); plot_error(stats_ewma.Z); axis tight; title('EWMA');
    drawnow
end

% -----------------------------------------------------------------------
% 2nd-level weighted (hierarchical) analysis
% -----------------------------------------------------------------------
[p,tm,Zcor,hewma_ste,hewma_mean,tvals,sb,stats_hewma] = hewma2(stats_ewma.Z, stats_ewma.var, ...
    lam, doplot, dodetrend, grpcontrast, base_timepts);

% -----------------------------------------------------------------------
% Mixture model to get classes
% -----------------------------------------------------------------------
[ind,ind2,stats_mixture_model] = Gaussian_mix(hewma_mean,50,base_timepts, 0, doplot);
ind = logical(ind); 
ind2 = logical(ind2);

STATS = struct('hewma_mean', hewma_mean, 'hewma_ste', hewma_ste, 'baseline_index', ind, 'active_index', ind2, ...
    'stats_ewma', stats_ewma, 'stats_hewma', stats_hewma, 'stats_mixture_model', stats_mixture_model);

% -----------------------------------------------------------------------
% Another final plot
% -----------------------------------------------------------------------
if doplot
    create_figure('Timeseries', 3, 1, 1);
    subplot(3, 1, 3);
    create_figure('Group Timeseries');

    xvals = (1:length(hewma_mean)) ./ samprate;

    plot(xvals, hewma_mean, 'k', 'LineWidth', 2); axis tight
    xlabel('Time');
    ylabel('Mean EWMA stat');
    
    fill_around_line(hewma_mean, hewma_ste, 'k', xvals);
    plot_horizontal_line(0, 'k')
    
    uncor_p_thresh = .05;

    tthr = tinv(1 - (uncor_p_thresh / 2), stats_hewma.df);

    sig_vector = abs(stats_hewma.tvals) > tthr;
    sigy = NaN * zeros(size(hewma_mean));
    sigy(sig_vector) = hewma_mean(sig_vector);
    plot(xvals, 0 * sigy, 'gx-', 'LineWidth', 3);

    axis tight

    % add mixture model
    %ylim = [min(hewma_mean) - hewma_ste max(hewma_mean) + hewma_ste];
    
    [y_baseclass, y_activeclass] = deal(NaN * zeros(size(hewma_mean)));
    y_baseclass(ind) = hewma_mean(ind);
    y_activeclass(ind2) = hewma_mean(ind2);
    
    plot(xvals, y_activeclass, 'r', 'LineWidth', 2);
end

end
