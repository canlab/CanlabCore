function obj = bootstrap(obj, X, Y, varargin)
% bootstrap  Bootstrap the model to get weight-map stability and FDR thresholds.
%
% :Usage:
% ::
%     pm = bootstrap(pm, X, Y);
%     pm = bootstrap(pm, X, Y, 'nboot', 5000);
%     pm = bootstrap(pm, X, Y, 'groups', subject_id);
%
% Pipeline:
%   1. Pre-fit bad-data check; remove bad rows/cols.
%   2. If obj is not yet fitted, run fit() on the full data so the
%      canonical weights.w / ml_model are populated.
%   3. For each of nboot samples:
%        - resample observations with replacement (or whole groups,
%          if 'groups' is provided)
%        - fit a fresh clone(obj) on the resample
%        - record the weight vector
%   4. Compute per-feature mean, SE, Z, two-tailed p, FDR threshold,
%      FDR-significant mask, and FDR-thresholded weights.
%
% After:
%   obj.weights.boot_w        [p x nboot] bootstrap weight samples
%   obj.weights.boot_w_mean   [p x 1]    mean across bootstraps
%   obj.weights.boot_w_ste    [p x 1]    SD across bootstraps
%   obj.weights.z             [p x 1]    boot_w_mean / boot_w_ste
%   obj.weights.p             [p x 1]    two-tailed normal-approximation p
%   obj.weights.fdr_thr       scalar     FDR(p, .05) threshold
%   obj.weights.fdr_sig       [p x 1]    logical, p <= fdr_thr
%   obj.weights.thresh_fdr    [p x 1]    obj.weights.w masked by fdr_sig

    pi = inputParser; pi.KeepUnmatched = true;
    addParameter(pi, 'nboot',   []);
    addParameter(pi, 'groups',  []);
    addParameter(pi, 'verbose', true);
    parse(pi, varargin{:});

    nboot   = pi.Results.nboot;
    if isempty(nboot)
        nboot = obj.nboot;
        if isempty(nboot), nboot = 1000; end
    end
    groups  = pi.Results.groups;
    verbose = pi.Results.verbose;

    % 1. Bad data
    [oc, of] = predictive_model.detect_bad_data(X, Y);
    if any(oc) || any(of)
        X(oc, :) = []; Y(oc) = []; X(:, of) = [];
        if ~isempty(groups), groups(oc) = []; end
    end
    Y = Y(:);
    n = numel(Y);
    p_feat = size(X, 2);

    if ~isempty(obj.random_state), rng(obj.random_state); end

    % 2. Ensure canonical fit.
    if ~obj.is_fitted || isempty(obj.weights) || ~isfield(obj.weights, 'w') ...
            || isempty(obj.weights.w)
        obj = fit(obj, X, Y, 'id', groups);
    end

    % 3. Bootstrap loop
    boot_w = nan(p_feat, nboot);

    if verbose, fprintf('bootstrap: %d samples', nboot); end
    every = max(1, ceil(nboot / 20));
    for b = 1:nboot
        if verbose && mod(b, every) == 0, fprintf('.'); end

        if isempty(groups)
            % observation-level resampling with replacement
            idx = randi(n, [n, 1]);
        else
            % grouped resampling: sample whole subjects
            uniq = unique(groups);
            new_g = uniq(randi(numel(uniq), [numel(uniq), 1]));
            idx = cell(numel(new_g), 1);
            for g = 1:numel(new_g)
                idx{g} = find(groups == new_g(g));
            end
            idx = vertcat(idx{:});
        end

        try
            m = clone(obj);
            m = fit(m, X(idx, :), Y(idx));
            w_b = m.weights.w;
            if numel(w_b) == p_feat
                boot_w(:, b) = w_b(:);
            end
        catch
            % skip failed bootstrap (e.g. degenerate resample)
        end
    end
    if verbose, fprintf(' done\n'); end

    % 4. Inference
    boot_w_mean = nanmean(boot_w, 2); %#ok<NANMEAN>
    boot_w_ste  = nanstd(boot_w, 0, 2); %#ok<NANSTD>
    % Treat exactly-zero SE as missing (feature was constant across
    % bootstraps) — gives z = 0, p = 1, no significance.
    boot_w_ste(boot_w_ste == 0) = Inf;
    z     = boot_w_mean ./ boot_w_ste;
    pvals = 2 * (1 - normcdf(abs(z)));
    % Clamp p-values away from 0 so FDR.m doesn't treat exact zeros
    % as "ineligible". eps is ~2.2e-16, smaller than any meaningful p.
    pvals = max(pvals, eps);

    if exist('FDR', 'file') == 2
        fdr_thr = FDR(pvals, 0.05);
    else
        fdr_thr = [];
    end
    if isempty(fdr_thr), fdr_thr = -Inf; end
    fdr_sig = pvals <= fdr_thr;

    thresh = obj.weights.w;
    if numel(thresh) == p_feat
        thresh(~fdr_sig) = 0;
    end

    obj.weights.boot_w      = boot_w;
    obj.weights.boot_w_mean = boot_w_mean;
    obj.weights.boot_w_ste  = boot_w_ste;
    obj.weights.z           = z;
    obj.weights.p           = pvals;
    obj.weights.fdr_thr     = fdr_thr;
    obj.weights.fdr_sig     = fdr_sig;
    obj.weights.thresh_fdr  = thresh;

    obj.history{end+1, 1} = sprintf('bootstrap: %d samples, %d FDR-sig features (thr=%.4g)', ...
        nboot, sum(fdr_sig), fdr_thr);
end
