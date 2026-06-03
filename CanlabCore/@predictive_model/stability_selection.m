function obj = stability_selection(obj, X, Y, varargin)
% stability_selection  Per-bootstrap top-k feature selection; report frequency.
%
% Stability selection runs a bootstrap and asks, on each resample,
% "which features got the largest |weights|?". Features that
% reliably appear in the top-k across bootstraps are "stable" —
% the classifier picks them regardless of which subjects it sees.
%
% This is the appropriate inference for high-dimensional regularised
% linear models (e.g. linear_svm under xval_SVM(highdimensional)),
% where the ordinary bootstrap z/p collapses because the
% regularisation pins the weights to nearly the same solution every
% time. Stability selection asks a different question — not "is the
% weight non-zero" but "is the weight in the top-k consistently" —
% and gives sharper inference in that regime.
%
% Ref: Meinshausen & Buhlmann, "Stability Selection", J. R. Stat. Soc. B (2010).
%
% :Usage:
% ::
%     pm = stability_selection(pm, X, Y, 'nboot', 100, 'k', 200);
%     pm = stability_selection(pm, X, Y, 'nboot', 200, 'k', 200, ...
%                              'threshold', 0.6, 'groups', subject_id);
%
% :Optional Inputs (name/value):
%   'nboot'      number of bootstrap resamples (default 100, or obj.nboot)
%   'k'          top-k features by |w| per bootstrap (default = 10% of
%                features, with a floor of 10)
%   'threshold'  selection-frequency threshold for "stable" (default 0.6)
%   'groups'     grouping vector for subject-grouped resampling
%   'verbose'    default true
%
% :Inputs:
%
%   **obj:**  a @predictive_model with a linear algorithm.
%   **X:**    [n x p] predictor matrix.
%   **Y:**    [n x 1] outcome vector.
%
% :Outputs:
%
%   **obj:**
%        the @predictive_model with obj.diagnostics.stability_selection
%        populated: .nboot, .k, .threshold, .selection_count [p x 1],
%        .selection_freq [p x 1] in [0,1], .stable (logical mask),
%        .n_stable, .valid_boots. Map .selection_freq to voxel space by
%        stashing it as a weight vector and calling weight_image.
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm');
%     X = dat.dat'; Y = dat.Y; id = dat.metadata_table.subj_id;
%     pm = predictive_model('algorithm','linear_svm','task','classification');
%     pm = stability_selection(pm, X, Y, 'nboot', 100, 'k', 2000, ...
%                              'threshold', 0.6, 'groups', id);
%     pm.diagnostics.stability_selection.n_stable
%
% :See also:
%   bootstrap, select_features, weight_image, grid_search

    pi = inputParser; pi.KeepUnmatched = true;
    addParameter(pi, 'nboot',     []);
    addParameter(pi, 'k',         []);
    addParameter(pi, 'threshold', 0.6);
    addParameter(pi, 'groups',    []);
    addParameter(pi, 'verbose',   true);
    parse(pi, varargin{:});

    nboot     = pi.Results.nboot;
    if isempty(nboot)
        nboot = obj.nboot;
        if isempty(nboot), nboot = 100; end
    end
    k         = pi.Results.k;
    threshold = pi.Results.threshold;
    groups    = pi.Results.groups;
    verbose   = pi.Results.verbose;

    % 1. Bad-data
    [oc, of] = predictive_model.detect_bad_data(X, Y);
    if any(oc) || any(of)
        X(oc, :) = []; Y(oc) = []; X(:, of) = [];
        if ~isempty(groups), groups(oc) = []; end
    end
    Y      = Y(:);
    n      = numel(Y);
    p_feat = size(X, 2);

    if isempty(k)
        k = max(10, round(p_feat * 0.1));
    end
    k = min(k, p_feat);

    if ~isempty(obj.random_state), rng(obj.random_state); end

    % 2. Ensure fitted for canonical weights / algorithm.
    if ~obj.is_fitted || isempty(obj.weights) || ~isfield(obj.weights, 'w') ...
            || isempty(obj.weights.w)
        obj = fit(obj, X, Y, 'id', groups);
    end

    % 3. Bootstrap loop, counting top-k membership per feature
    selection_count = zeros(p_feat, 1);
    valid_boots     = 0;

    if verbose, fprintf('stability_selection: %d boots, top-k=%d', nboot, k); end
    every = max(1, ceil(nboot / 20));

    for b = 1:nboot
        if verbose && mod(b, every) == 0, fprintf('.'); end

        if isempty(groups)
            idx = randi(n, [n, 1]);
        else
            uniq    = unique(groups);
            new_g   = uniq(randi(numel(uniq), [numel(uniq), 1]));
            idx_c   = cell(numel(new_g), 1);
            for g = 1:numel(new_g)
                idx_c{g} = find(groups == new_g(g));
            end
            idx = vertcat(idx_c{:});
        end

        try
            m = clone(obj);
            m = fit(m, X(idx, :), Y(idx));
            w = m.weights.w;
            if numel(w) == p_feat && ~all(isnan(w))
                [~, ord] = sort(abs(w), 'descend');
                top_k = ord(1:k);
                selection_count(top_k) = selection_count(top_k) + 1;
                valid_boots = valid_boots + 1;
            end
        catch
            % skip degenerate resample
        end
    end
    if verbose, fprintf(' done\n'); end

    if valid_boots == 0
        warning('predictive_model:stability_selection:NoValidBoots', ...
            'All %d bootstrap fits failed; stability_selection results empty.', nboot);
        selection_freq = nan(p_feat, 1);
    else
        selection_freq = selection_count / valid_boots;
    end
    stable = selection_freq >= threshold;

    obj.diagnostics.stability_selection = struct( ...
        'nboot',           nboot, ...
        'k',               k, ...
        'threshold',       threshold, ...
        'selection_count', selection_count, ...
        'selection_freq',  selection_freq, ...
        'stable',          stable, ...
        'n_stable',        sum(stable), ...
        'valid_boots',     valid_boots);

    obj.history{end+1, 1} = sprintf( ...
        'stability_selection: %d valid boots, top-k=%d, thr=%.2f -> %d stable features', ...
        valid_boots, k, threshold, sum(stable));
end
