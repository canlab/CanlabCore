function obj = select_features(obj, X, Y, varargin)
% select_features  Univariate feature selection.
%
% Computes a per-feature p-value against Y (two-sample t-test for binary
% classification, one-way ANOVA for multiclass, Pearson correlation
% p-value for regression). Selects features by p-value threshold or
% by top-k.
%
% :Usage:
% ::
%     pm = select_features(pm, X, Y);                 % default: p < .05
%     pm = select_features(pm, X, Y, 'pthresh', .01);
%     pm = select_features(pm, X, Y, 'k', 100);       % top 100 by p
%     pm = select_features(pm, X, Y, 'method', 'top_k_correlation', 'k', 50);
%
% The selection is stored in obj.diagnostics.feature_selection and is
% also unioned into obj.omitted_features (so a subsequent predict(pm,
% X_new) on full-width X_new will automatically drop the non-selected
% columns). If you fit() after select_features, fit() will detect bad
% data afresh and overwrite omitted_features — call select_features
% AFTER fit for masks to persist.
%
% IMPORTANT: when used in conjunction with crossval, the selection is
% computed on the WHOLE dataset, which leaks held-out information into
% feature choice. For honest CV, run select_features per fold (manually
% or via a future select_features-inside-Pipeline composer).
%
% :Optional Inputs (name/value):
%   'method'  'univariate' (default; t-test / ANOVA / correlation) or
%             'top_k_correlation' (rank by |corr(x_j, Y)|)
%   'pthresh' p-value threshold (default 0.05; ignored if 'k' set)
%   'k'       keep top-k features by score (overrides 'pthresh')
%   'verbose' default true
%
% After:
%   obj.omitted_features                   (unioned with selection mask)
%   obj.diagnostics.feature_selection      struct:
%     .method      char
%     .pvalues     [p x 1] (univariate only)
%     .scores      [p x 1] (top_k_correlation: |corr|)
%     .selected    [p x 1] logical
%     .n_selected, .n_input

    pi = inputParser; pi.KeepUnmatched = true;
    addParameter(pi, 'method',  'univariate');
    addParameter(pi, 'pthresh', 0.05);
    addParameter(pi, 'k',       []);
    addParameter(pi, 'verbose', true);
    parse(pi, varargin{:});

    method  = lower(pi.Results.method);
    pthresh = pi.Results.pthresh;
    k       = pi.Results.k;
    verbose = pi.Results.verbose;

    n_features_in = size(X, 2);
    Y = Y(:);

    if isempty(obj.task)
        if numel(unique(Y(~isnan(Y)))) <= 2
            obj.task = 'classification';
        else
            obj.task = 'regression';
        end
    end

    pvals  = [];
    scores = [];

    switch method
        case 'univariate'
            pvals = compute_univariate_pvals(X, Y, obj.task);
            if ~isempty(k)
                [~, idx] = sort(pvals);
                keep = false(n_features_in, 1);
                keep(idx(1:min(k, n_features_in))) = true;
            else
                keep = pvals(:) < pthresh;
            end

        case 'top_k_correlation'
            if isempty(k)
                error('predictive_model:select_features:NoK', ...
                    'top_k_correlation requires ''k''.');
            end
            scores = zeros(n_features_in, 1);
            for j = 1:n_features_in
                scores(j) = abs(corr(X(:, j), Y, 'rows', 'complete'));
            end
            [~, idx] = sort(scores, 'descend');
            keep = false(n_features_in, 1);
            keep(idx(1:min(k, n_features_in))) = true;

        otherwise
            error('predictive_model:select_features:UnknownMethod', ...
                'Unknown method ''%s''. Options: ''univariate'', ''top_k_correlation''.', method);
    end

    % Union into omitted_features (drop non-selected on top of any existing mask).
    if isempty(obj.omitted_features) || ~islogical(obj.omitted_features) ...
            || numel(obj.omitted_features) ~= n_features_in
        obj.omitted_features = false(n_features_in, 1);
    end
    obj.omitted_features = obj.omitted_features | ~keep(:);

    obj.diagnostics.feature_selection = struct( ...
        'method',     method, ...
        'pvalues',    pvals, ...
        'scores',     scores, ...
        'selected',   keep, ...
        'n_selected', sum(keep), ...
        'n_input',    n_features_in);

    obj.history{end+1, 1} = sprintf( ...
        'select_features (%s): kept %d / %d features', ...
        method, sum(keep), n_features_in);

    if verbose
        fprintf('select_features (%s): kept %d / %d features\n', ...
            method, sum(keep), n_features_in);
    end
end


% --- helpers ---------------------------------------------------------------
function pvals = compute_univariate_pvals(X, Y, task)
    [~, p] = size(X);
    pvals = nan(p, 1);

    if strcmp(task, 'classification')
        classes = unique(Y(~isnan(Y)));
        if numel(classes) == 2
            msk1 = Y == classes(1);
            msk2 = Y == classes(2);
            for j = 1:p
                [~, pvals(j)] = ttest2(X(msk1, j), X(msk2, j));
            end
        else
            % multiclass: one-way ANOVA per feature
            for j = 1:p
                try
                    pvals(j) = anova1(X(:, j), Y, 'off');
                catch
                    pvals(j) = NaN;
                end
            end
        end
    else
        % regression: Pearson correlation p-value
        for j = 1:p
            [~, pvals(j)] = corr(X(:, j), Y, 'rows', 'complete');
        end
    end
end
