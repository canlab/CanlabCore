function stats = xcorr_xy_multisubject(X, Y)
% This function will cross-correlate two variables, X and Y for each of N
% subjects.  Correlation and latency values will be saved.
%
% Second-level tests are done across the N subjects on each of the
% correlation and latency values.
%
% :Inputs:
%
%   **X:**
%        must be an observations x N matrix
%
%        NOTE: observations are assumed to be timeseries values!
%
%        this DOES matter because an AR(2) model is used...see below
%
%   **Y:**
%        must be the same.
%
%        cross-correlations will be computed for pairs of columns, separately for 
%        each successive column of X / Y
%
% This function ues shift_correl.m, with a two-pass procedure.
%
% The first pass provides initial estimates, including an estimate of the latency standard
% deviation, which is used as an Empirical Bayes prior.  The second pass is used to
% apply the EBayes priors and estimate latencies and cross-correlation values.
% An ar(2) model is used to estimate the DF.
%
% The estimates will be biased towards zero in the second pass, but they
% will likely be lower variance.
%
% stats.latency and stats.association return group-level stats on these
% computed using a sign pemutation test (to avoid normality assumption).
% the "association" test is a test of the relationship, but is different
% than a standard multilevel model because it uses correlation values rather
% than slope values, and the sign permutation test.

detailplots = 1;

% -------------------------------------------------------------------------
% FIRST PASS
% -------------------------------------------------------------------------

for i = 1:size(X, 2)
        if ~any(isnan(Y(:, i)))
            [shiftvals(i, :), corrvals(i, :), bestshift(i, 1), bestcorr(i, 1)] = ...
                shift_correl(X(:, i), Y(:, i), 'max_shift', 4, 'shift_step', .2, 'use_rob_stats', 0, ...
                'return_betas', 0, 'display_plots', 0,'optimize');
        else
            shiftvals(i, :) = NaN;
            corrvals(i, :) = NaN;
            [bestshift(i, 1), bestcorr(i, 1)] = deal(NaN);
        end
end

    
% -------------------------------------------------------------------------
% SECOND PASS
% -------------------------------------------------------------------------

    ebayes_priors = [0 nanstd(bestshift)];
    
    for i = 1:size(X, 2)
        if ~any(isnan(Y(:, i)))
            [shiftvals(i, :), corrvals(i, :), bestshift(i, 1), bestcorr(i, 1), aout, bout] = ...
                shift_correl(X(:, i), Y(:, i), 'max_shift', 4, 'shift_step', .2, 'use_rob_stats', 0, ...
                'return_betas', 0, 'display_plots', detailplots, 'optimize', 'priors', ebayes_priors, 'ar2');
            close
        else
            shiftvals(i, :) = NaN;
            corrvals(i, :) = NaN;
            [bestshift(i, 1), bestcorr(i, 1)] = deal(NaN);
        end
    end
    
    
% -------------------------------------------------------------------------
% STATS
% -------------------------------------------------------------------------
X2 = ones(size(bestshift));
stats.latency = glmfit_general(bestshift, X2, 'signperm', 'nresample', 10000);

stats.association = glmfit_general(bestcorr, X2, 'signperm', 'nresample', 10000);

stats.latency.individuals = bestshift;
stats.association.individuals = bestcorr;

% -------------------------------------------------------------------------
% PLOT
% -------------------------------------------------------------------------

    create_figure('Shift', 1, 2); hist(bestshift, 30); plot_vertical_line(mean(bestshift));
    subplot(1, 2, 2); plot_error(mean(shiftvals)', corrvals);
    hold on
    plot_vertical_line(0);
    set(gca, 'FontSize', 24);
    xlabel('Time shift');
    ylabel('Cross-correlation');
    axis auto
    axis tight
    
end
