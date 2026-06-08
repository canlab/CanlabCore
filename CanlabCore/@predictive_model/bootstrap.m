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
%   obj.weights.z             [p x 1]    boot_w_mean / boot_w_ste (SE floored
%                                        so it can't produce 1e15 values
%                                        from numerical jitter)
%   obj.weights.p             [p x 1]    continuity-corrected empirical
%                                        two-tailed bootstrap p:
%                                          2 * (min(n_pos, n_neg) + 1)
%                                              / (n_valid + 1)
%                                        bounded in [2/(nboot+1), 1].
%                                        This is more robust than the
%                                        Wald-style z->p mapping when
%                                        bootstrap distributions are
%                                        non-Gaussian or numerically
%                                        constant (the common case for
%                                        regularised linear SVM).
%   obj.weights.p_wald        [p x 1]    Wald-style 2*(1 - normcdf(|z|))
%                                        kept for sanity-checking / continuity
%                                        with legacy CANlab code.
%   obj.weights.fdr_thr       scalar     FDR(p, .05) threshold on the
%                                        empirical p above
%   obj.weights.fdr_sig       [p x 1]    logical, p <= fdr_thr
%   obj.weights.thresh_fdr    [p x 1]    obj.weights.w masked by fdr_sig
%
% NOTE on regularised linear SVM (e.g. fitclinear under
% xval_SVM('highdimensional', true)): bootstrap weights are often
% nearly numerically identical across resamples because L2 regularisation
% drives the optimiser to almost the same solution regardless of the
% specific sample. In that regime the SE is dominated by numerical
% jitter, the Wald z explodes, and bootstrap-based FDR will flag almost
% every voxel. The empirical p above is more honest — it bottoms out
% at 2/(nboot+1), which is the smallest p any bootstrap procedure can
% give you given the resolution of the resampling. If you want
% finer-grained voxel-wise inference, fit a less-aggressively-regularised
% model (smaller Lambda) or use permutation testing on a univariate
% statistic.
%
% :Inputs:
%
%   **obj:**
%        a @predictive_model with obj.algorithm set (need not be pre-fit).
%
%   **X:**
%        [n x p] predictor matrix.
%
%   **Y:**
%        [n x 1] outcome vector.
%
% :Optional Inputs (name/value):
%   'nboot'    number of bootstrap resamples (default obj.nboot or 1000)
%   'groups'   grouping vector; resample whole groups (clusters) instead
%              of individual observations (e.g. subject id)
%
% :Outputs:
%
%   **obj:**
%        the @predictive_model with bootstrap weight statistics in
%        obj.weights (.boot_w, .boot_w_mean, .boot_w_ste, .z, .p, .p_wald,
%        .fdr_thr, .fdr_sig, .thresh_fdr). Map to voxel space with
%        weight_map_object(obj, source).
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm');
%     X = dat.dat'; Y = dat.Y; id = dat.metadata_table.subj_id;
%     pm = predictive_model('algorithm','svm','task','classification');
%     pm = crossval(pm, X, Y, 'groups', id);
%     pm = bootstrap(pm, X, Y, 'nboot', 1000, 'groups', id);
%     sum(pm.weights.fdr_sig)          % # FDR-significant voxels
%     [~, si] = weight_map_object(pm, dat);   % statistic_image with .p / .sig
%
% :See also:
%   weight_map_object, crossval, permutation_test, stability_selection

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

    % Normalize non-numeric groups (cell of subject-id strings, etc.).
    if ~isempty(groups) && ~(isnumeric(groups) || islogical(groups))
        [~, ~, groups] = unique(groups(:), 'stable');
    end

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

    % Z statistic. Floor the SE before division so numerical jitter
    % (e.g. fitclinear converging to the same w +- 1e-20 each
    % bootstrap) can't produce z = 1e15. Floor is the larger of
    % eps('single') in absolute terms and 1e-6 * |mean| relative.
    abs_mean   = abs(boot_w_mean);
    se_floor   = max(eps('single'), 1e-6 * abs_mean);
    ste_safe   = max(boot_w_ste, se_floor);
    z          = boot_w_mean ./ ste_safe;
    p_wald     = 2 * (1 - normcdf(abs(z)));
    p_wald     = max(p_wald, eps);

    % Empirical bootstrap two-tailed p-value with continuity correction.
    % For each feature, ask: what fraction of bootstrap samples landed
    % on the "wrong" side of zero? Doubled for two-tailed. Naturally
    % bounded in [2/(nboot+1), 1] so it can't collapse to numerical
    % zero, and it doesn't depend on normality.
    n_valid    = sum(~isnan(boot_w), 2);
    n_pos      = sum(boot_w > 0, 2);
    n_neg      = sum(boot_w < 0, 2);
    pvals_emp  = 2 * (min(n_pos, n_neg) + 1) ./ (n_valid + 1);
    pvals_emp  = min(pvals_emp, 1);

    % FDR threshold runs on the empirical p (default; more honest than
    % p_wald when bootstrap is numerically pinned).
    if exist('FDR', 'file') == 2
        fdr_thr = FDR(pvals_emp, 0.05);
    else
        fdr_thr = [];
    end
    if isempty(fdr_thr), fdr_thr = -Inf; end
    fdr_sig = pvals_emp <= fdr_thr;

    thresh = obj.weights.w;
    if numel(thresh) == p_feat
        thresh(~fdr_sig) = 0;
    end

    obj.weights.boot_w      = boot_w;
    obj.weights.boot_w_mean = boot_w_mean;
    obj.weights.boot_w_ste  = boot_w_ste;
    obj.weights.z           = z;
    obj.weights.p           = pvals_emp;
    obj.weights.p_wald      = p_wald;
    obj.weights.fdr_thr     = fdr_thr;
    obj.weights.fdr_sig     = fdr_sig;
    obj.weights.thresh_fdr  = thresh;

    obj.history{end+1, 1} = sprintf( ...
        'bootstrap: %d samples, %d FDR-sig features (empirical p, thr=%.4g)', ...
        nboot, sum(fdr_sig), fdr_thr);
end
