function [shiftvals, corrvals, bestshift, bestcorr, aout, bout] = shift_correl(a, b, varargin)
% :Usage:
% ::
%
%     [shiftvals, corrvals, bestshift, bestcorr, aout, bout] = shift_correl(a, b, ['max_shift', max_shift], ['shift_step', shift_step], ['use_rob_stats', 0|1], ['return_betas', 0|1], ['display_plots', 0|1],['priors',[mu sig]])
%
% Shifts a backwards and forwards by max_shift elements (default 12)
% in each direction, computing the correlation
% between a and b at each shift value
%
% Thus, negative values mean a happens AFTER b
% Positive values means b happens later
%
% truncates the tails of a and b where necessary
% to ensure they're the same length.
%
% :Inputs:
%
%   **a, b:**
%        two vectors to correlate
%
%   **max_shift:**
%        max number of elements to shift (default 12)
%
%   **shift_step:**
%        size of incremental shifts (default .1)
%
%   **use_rob_stats:**
%        robust regression, IRLS (default 0)
%
%   **return_betas:**
%        return betas rather than corr coeffs (default 0)
%
%   **display_plots:**
%        display plot of shift vals and correlations/betas (default 0)
%
%        NOTE: robust r value of chosen solution may not be max,
%        as weighted_corrcoef now does not account for variance
%        inflation with low weights for some observations
%
%   **priors:**
%        Incorporate gaussian priors with mean priors(1) and std
%        priors(2)
%
% :Examples:
% ::
%
%    % first extract data for a subject into clusters:
%    clusters = roi_probe(spm_get(Inf), 'snpm0002_varsm_cov_clusters.mat');
%
%    % OR
%
%    clusters = mask2clusters('SnPMt_filtered.img', spm_get(Inf));
%
%    % Then do the shifting:
%
%    for i=1:length(clusters)
%        [shiftvals, corrvals] = shift_correl(xX.X(:, 1), clusters(i).timeseries);
%        title(['Cluster' num2str(i)]);
%        disp(['Cluster ' num2str(i) ' estimate shift by ' num2str(1.5 * shiftvals(corrvals == max(corrvals)))]);
%    end
%
%
%    % Simulation
%    trueval = 1;
%    n = 100;
%    noisevar = 0;
%    a = randn(n,1);
%    b = a + randn(n,1).*noisevar;
%    a=smooth_timeseries(a,n./5);
%    b=smooth_timeseries(b,n./5);
%    a = shift_signal(a,trueval);
%    b = b(1:length(a));
%    figure;plot(a);
%    hold on;
%    plot(b,'r')
%    [shiftvals, corrvals, bestshift, bestcorr, aout, bout] = ...
%    shift_correl(a, b, 'max_shift', 4, 'shift_step', .2, 'use_rob_stats', 1, 'return_betas', 0, 'display_plots', 1,'optimize'); bestshift
%
%    r = .7; shiftstep = 1;
%    optstr = 'optimize'; %    or 'noopt'
%    tic
%    for i = 1:10
%        ab = mvnrnd([0 0],[1 r; r 1],n);
%        a = ab(:,1);
%        b = ab(:,2);
%        a=smooth_timeseries(a,n./5);
%        b=smooth_timeseries(b,n./5);
%        a = shift_signal(a,trueval);
%        b = b(1:length(a));
%        [shiftvals, corrvals, bestshift(i), bestcorr(i), aout, bout] = shift_correl(a, b, 'max_shift', 4, 'shift_step', shiftstep, 'use_rob_stats', 1,optstr);
%    end
%    toc

    displaying_plots = 0;
    using_rob_stats = 0;
    returning_betas = 0;
    shiftvals = [];
    corrvals = [];
    max_shift = 12;
    shift_step = 1;
    do_optimize = 0;
    priors = [];
ar_model = 0;

    for i=1:length(varargin)
        if ischar(varargin{i})
            switch(varargin{i})
                case 'max_shift'
                    max_shift = varargin{i+1};
                case 'shift_step'
                    shift_step = varargin{i+1};
                case {'using_rob_stats' 'use_rob_stats'}
                    warning('off','stats:statrobustfit:IterationLimit');
                    using_rob_stats = varargin{i+1};
                case {'returning_betas' 'return_betas'}
                    returning_betas = varargin{i+1};
                case {'displaying_plots' 'display_plots'}
                    displaying_plots = varargin{i+1};
                case 'optimize'
                    do_optimize = 1;
                case {'prior' 'priors' 'prior_prob'}
                    priors = varargin{i+1};
                    
                case {'ar', 'ar1', 'ar_model'}
                    ar_model = 1;
                    
                case {'ar2'}
                    ar_model = 2;
            end
        end
    end

    if(size(a,2)~=1), a = a'; end
    if(size(b,2)~=1), b = b'; end

    whnan = find(isnan(a) | isnan(b));
    a(whnan) = [];
    b(whnan) = [];

    shiftvals = -max_shift:shift_step:max_shift;
    len = length(shiftvals);
    corrvals = zeros(1,len);
    pvals = zeros(1,len);

    for shift = 1:len

        pvals(shift) = do_shift_corr(shiftvals(shift),a,b,[],using_rob_stats,returning_betas,priors, ar_model);

    end

    bestidx = find(pvals == min(pvals));
    bestshift = shiftvals(bestidx);

    % get aout and bout, based on bestshift
    truncate_output_vectors();

    bestcorr = compute_corrs_or_betas(aout,bout,using_rob_stats,returning_betas);
    %bestcorr = corrvals(bestidx);

    % refine search with nonlinear fit
    if do_optimize
        options = optimset('MaxFunEvals',100000,'Maxiter',10000000,'Display','off','LevenbergMarquardt','on');
        limits = [bestshift - .99 bestshift + .99];
        bestshift = fminsearch(@(shift) do_shift_corr(shift,a,b,limits,using_rob_stats,returning_betas,priors,ar_model),bestshift,options);

        % get aout and bout, based on bestshift
        truncate_output_vectors();

        % get best corr
        bestcorr = compute_corrs_or_betas(aout,bout,using_rob_stats,returning_betas);
    else
        truncate_output_vectors();
    end

    if displaying_plots
        % get all corrvals
        for shift = 1:len
            [shiftedvec,unshiftedvec] = shift_vec(a,b,shiftvals(shift));
            corrvals(shift) = compute_corrs_or_betas(shiftedvec,unshiftedvec,using_rob_stats,returning_betas);
        end

        display_plots();
    end





    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Nested functions
    %%%%%%%%%%%%%%%%%%%%%%%%%



    function truncate_output_vectors()
        if bestshift > 0
            bout = shift_signal(b, bestshift);
            aout = a(1:length(bout));
        elseif bestshift < 0
            aout = shift_signal(a, bestshift);
            bout = b(1:length(aout));
        else
            aout = a;
            bout = b;
        end
    end


    function display_plots()
        create_figure('shift_correl_plot', 1, 2);

        subplot(1, 2, 1);
        plot(shiftvals, corrvals, 'ko-', 'LineWidth', 2); hold on;
        plot(bestshift, bestcorr, 'rx', 'MarkerSize', 16, 'LineWidth', 2);

        xlabel('Shift (in elements)');
        if(returning_betas)
            ylabel('Betas');
        else
            ylabel('Correlation');
        end

        subplot(1, 2, 2);
        plot_correlation_samefig(aout, bout, [], 'k.', 0, using_rob_stats);
        xlabel('a');
        ylabel('b');
    end
    
end  % main function


% this function is used in fminsearch optimization
% returns p-values, which are minimized in the optimization
% alternatively, if return betas is specified, returns 1./abs(beta) which
% should also be minimized by fminsearch or comparable algorithm
function pval = do_shift_corr(shift,a,b,limits,using_rob_stats,returning_betas,priors,ar_model)

    % constrain to prespecified limits
    if ~isempty(limits), shift = max(shift,limits(1)); shift = min(shift,limits(2)); end

    %if abs(shift) > max_shift, shift = sign(shift) .* max_shift; end

    % compute
    [shiftedvec,unshiftedvec] = shift_vec(a,b,shift);
    compute_pvals_or_betas();  % returns log pval, or 1/abs(slope beta).  Lower values are better.

    % incorporate priors, if we have them
    % P(h | dat) = P(dat | h)P(dat) / P(h)
    % assume P(dat) = 1.  It's unknown.  P(h) is priorp
    if ~isempty(priors)
        priorp = log(max(normpdf(shift,priors(1),priors(2)),eps));
        pval = pval - priorp; % This is P(dat) / P(h), in log units. pval is actually log. adding logs multiplies p-values.
    end

    %correrr = 1 ./ abs(corrval);

    %inline
    % -----------------------------------
    function compute_pvals_or_betas()
        % Returns subtractand of p-values ( 1 - p ) or betas
        if all(shiftedvec - mean(shiftedvec) < eps) || all(unshiftedvec - mean(unshiftedvec) < eps)
            % no variance in a or b
            pval = NaN;
        elseif using_rob_stats && returning_betas
            betas = robustfit(shiftedvec, unshiftedvec);
            pval = 1./abs(betas(2));
            
        elseif using_rob_stats
            [betas, stats] =  robustfit(shiftedvec, unshiftedvec, 'bisquare'); %#ok
            pval = log(stats.p(2));
            %r = weighted_corrcoef([shiftedvec unshiftedvec], stats.w);
            %corrval = r(1, 2);
            
        elseif returning_betas
            betas = glmfit(shiftedvec, unshiftedvec);
            pval = 1./abs(betas(2));
            %corrval = betas(2);
            
        elseif ar_model  % ar(2)
            [betas, t, pval] = fit_gls(shiftedvec,unshiftedvec, [], ar_model);
            
        else
            [R,p] = corrcoef(shiftedvec, unshiftedvec);
            %corrval = R(1, 2);
            pval = log(p(1,2));
        end
    end
end



function [shiftedvec,unshiftedvec] = shift_vec(a,b,shift)
    if(shift > 0)
        shiftedvec = shift_signal(b, shift);
        unshiftedvec = a(1:length(shiftedvec));
    else
        shiftedvec = shift_signal(a, shift);
        unshiftedvec = b(1:length(shiftedvec));
    end
end


function val = compute_corrs_or_betas(shiftedvec,unshiftedvec,using_rob_stats,returning_betas)

    % remove for speed
    %if all(shiftedvec - mean(shiftedvec) < eps) || all(unshiftedvec - mean(unshiftedvec) < eps)
    % no variance in a or b
    % val = NaN;
    if using_rob_stats && returning_betas
        betas = robustfit(shiftedvec, unshiftedvec);
        val = betas(2);
    elseif using_rob_stats
        [betas, stats] =  robustfit(shiftedvec, unshiftedvec, 'bisquare'); %#ok
        %             if(~iscol(shiftedvec)), shiftedvec = shiftedvec'; end
        %             if(~iscol(unshiftedvec)), unshiftedvec = shiftedvec'; end
        r = weighted_corrcoef([shiftedvec unshiftedvec], stats.w);
        val = r(1, 2);
    elseif returning_betas
        betas = glmfit(shiftedvec, unshiftedvec);
        val = betas(2);
    else
        R = corrcoef(shiftedvec, unshiftedvec);
        val = R(1, 2);
    end
end
