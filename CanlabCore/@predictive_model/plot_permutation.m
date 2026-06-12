function varargout = plot_permutation(obj, varargin)
% plot_permutation  Histogram of the permutation null distribution.
%
% Draws the null distribution of the cross-validated score built by
% permutation_test, with the observed score marked by a solid line and the
% alpha-level critical value (the "significance threshold" the observed score
% must beat) marked by a dashed line. The shaded tail past the critical value
% is the rejection region.
%
% :Usage:
% ::
%     plot_permutation(pm);
%     plot_permutation(pm, 'alpha', 0.01);
%     h = plot_permutation(pm, 'nbins', 40);
%
% :Inputs:
%
%   **obj:**
%        a @predictive_model on which permutation_test has been run (so
%        obj.permutation_results is populated).
%
% :Optional Inputs (name/value):
%
%   **'alpha':**
%        significance level for the critical-value line (default 0.05).
%
%   **'nbins':**
%        number of histogram bins (default 30).
%
%   **'noplot':**
%        compute/return the summary struct without drawing.
%
% :Outputs:
%
%   **h:**
%        (optional) the axes handle (or, with 'noplot', a struct with the
%        observed score, p-value, and critical value).
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm', 'noverbose');
%     X = dat.dat'; Y = dat.Y; id = dat.metadata_table.subj_id;
%     pm = predictive_model('algorithm','svm','task','classification');
%     pm = crossval(pm, X, Y, 'groups', id);
%     pm = permutation_test(pm, X, Y, 'nperm', 1000, 'groups', id);
%     plot_permutation(pm);
%
% :See also:
%   permutation_test, plot, summary

    p = inputParser; p.KeepUnmatched = true;
    addParameter(p, 'alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
    addParameter(p, 'nbins', 30, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'noplot', false, @(x) islogical(x) || isnumeric(x));
    % accept bare 'noplot' flag
    raw = varargin;
    is_np = cellfun(@(a) (ischar(a) || isstring(a)) && strcmpi(a, 'noplot'), raw);
    do_noplot_flag = any(is_np);
    raw(is_np) = [];
    parse(p, raw{:});
    alpha = p.Results.alpha;
    nbins = p.Results.nbins;
    do_plot = ~(logical(p.Results.noplot) || do_noplot_flag);

    pr = obj.permutation_results;
    if ~isstruct(pr) || ~isfield(pr, 'null_scores') || isempty(pr.null_scores)
        error('predictive_model:plot_permutation:NoResults', ...
            'No permutation results. Run permutation_test(pm, X, Y) first.');
    end

    null_scores = pr.null_scores(:);
    null_scores = null_scores(~isnan(null_scores));
    obs = pr.observed;
    pval = pr.p_value;

    % Direction: greater-is-better scorers reject in the upper tail.
    gib = true;
    try
        sc = cv_scorer.make(pr.scorer);
        gib = logical(sc.greater_is_better);
    catch
    end

    if gib
        crit = quantile(null_scores, 1 - alpha);
    else
        crit = quantile(null_scores, alpha);
    end

    summary_struct = struct('observed', obs, 'p_value', pval, ...
        'critical_value', crit, 'alpha', alpha, 'scorer', pr.scorer, ...
        'permutation', pr.permutation, 'nperm', numel(null_scores));

    if do_plot
        ax = gca; cla(ax); hold(ax, 'on');
        histogram(ax, null_scores, nbins, 'FaceColor', [.6 .6 .6], ...
            'EdgeColor', 'none', 'FaceAlpha', .9);
        yl = ylim(ax);

        % critical value (dashed) and observed (solid).
        plot(ax, [crit crit], yl, '--', 'Color', [.85 .33 .10], 'LineWidth', 2);
        plot(ax, [obs obs],   yl, '-',  'Color', [0 0 0],       'LineWidth', 2.5);

        xlabel(ax, sprintf('Null %s', strrep(pr.scorer, '_', ' ')));
        ylabel(ax, 'Permutations');
        title(ax, sprintf('Permutation null (%s): observed = %.3g, p = %.4g', ...
            strrep(pr.permutation, '_', ' '), obs, pval));
        legend(ax, {'null', sprintf('p = %.2g threshold', alpha), 'observed'}, ...
            'Location', 'best'); legend(ax, 'boxoff');
        hold(ax, 'off');
        set(gcf, 'Color', 'w');
    end

    if nargout > 0
        if do_plot, varargout{1} = gca; else, varargout{1} = summary_struct; end
    end
end
