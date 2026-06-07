function results = compare_predict_legacy_vs_newapi(dat, varargin)
% compare_predict_legacy_vs_newapi  Legacy vs 'newapi' fmri_data.predict.
%
% Re-runnable diagnostic: runs fmri_data.predict the legacy way and again
% with the 'newapi' keyword, **holding the cross-validation fold
% assignment constant** (an explicit integer fold vector is passed to both
% so any difference is purely the algorithm, not the partition), and
% reports whether the cross-validated predictions and error match.
%
% :Usage:
% ::
%     results = compare_predict_legacy_vs_newapi();              % emotionreg
%     results = compare_predict_legacy_vs_newapi(my_fmri_data);
%     results = compare_predict_legacy_vs_newapi(dat, 'nfolds', 5);
%     results = compare_predict_legacy_vs_newapi(dat, 'tol', 1e-6);
%
% :Inputs:
%
%   **dat:**
%        an fmri_data object. If omitted, load_image_set('emotionreg').
%        Its .Y is overwritten per-algorithm: a balanced binary +/-1 Y for
%        classification algorithms, a continuous Y for regression
%        algorithms (so each algorithm is tested on an appropriate outcome
%        with the SAME images and SAME folds).
%
% :Optional Inputs (name/value):
%
%   **'nfolds':** scalar k for the fixed fold vector (default 5).
%   **'tol':**    max-abs-difference tolerance to call outputs "identical"
%                 (default 1e-6).
%   **'algorithms':** cell of algorithm specs to test. Each entry is
%                 either an algorithm_name string, or a 2-cell
%                 {algorithm_name, extra_opts_cell} where extra_opts_cell
%                 are extra name/value or flag args forwarded to BOTH the
%                 legacy and newapi predict() calls (e.g. {'lasso_num',5}).
%                 Default:
%                   {'cv_svm','cv_svr','cv_pcr','cv_lassopcr', ...
%                    {'cv_lassopcr',{'lasso_num',5}}, ...
%                    {'cv_lassopcr',{'estimateparam'}}}.
%   **'seed':**   rng seed for the synthetic outcomes / fold vector (default 1).
%
% :Outputs:
%
%   **results:** a table, one row per algorithm, with:
%        algorithm, task, max_abs_yfit_diff, corr_yfit, cverr_legacy,
%        cverr_newapi, cverr_absdiff, identical (logical at tol).
%
% Interpretation (with folds held constant):
%   - 'cv_svm'      -> 'svm' (fitcsvm linear): EXACT match.
%   - 'cv_pcr'      -> 'pcr' (PCA + OLS, faithful to legacy cv_pcr): EXACT.
%   - 'cv_lassopcr' (default)        -> 'pcr'     : EXACT (full model OLS-refit == PCR).
%   - 'cv_lassopcr' {'lasso_num',k}  -> 'lassopcr': EXACT (same LASSO path +
%                     relaxed-OLS refit on the non-zero components).
%   - 'cv_lassopcr' {'estimateparam'}-> 'lassopcr': runs nested-CV lambda
%                     selection. The new path uses a clean inner CV, so it
%                     can DIFFER from the legacy estimateparam value (the
%                     legacy nested CV reuses the outer fold vector); the new
%                     behavior is the corrected one. Not expected identical.
%   - 'cv_svr'      -> 'svr' (fitrsvm linear): NEAR-identical (corr ~0.9997),
%                     not exact, because the legacy Spider SVR solver and
%                     fitrsvm are different implementations. This is a solver
%                     difference, not a fold bug.
%
% :Examples:
% ::
%     r = compare_predict_legacy_vs_newapi();
%     disp(r);
%
% :See also:
%   fmri_data.predict, predictive_model, predictive_model_unit_test

    % Allow calling with only name/value pairs (no dat): shift if the first
    % positional is actually an option name.
    if nargin >= 1 && (ischar(dat) || isstring(dat))
        varargin = [{dat}, varargin];
        dat = [];
    end

    p = inputParser;
    addParameter(p, 'nfolds', 5);
    addParameter(p, 'tol', 1e-6);
    addParameter(p, 'algorithms', {'cv_svm','cv_svr','cv_pcr','cv_lassopcr', ...
        {'cv_lassopcr', {'lasso_num', 5}}, {'cv_lassopcr', {'estimateparam'}}});
    addParameter(p, 'seed', 1);
    parse(p, varargin{:});
    nfolds = p.Results.nfolds;
    tol    = p.Results.tol;
    algos  = p.Results.algorithms;
    seed   = p.Results.seed;

    if nargin < 1 || isempty(dat)
        dat = load_image_set('emotionreg', 'noverbose');
    end
    dat.dat = double(dat.dat);
    n = size(dat.dat, 2);

    rng(seed);

    % --- Fixed fold vector (held constant across legacy + newapi) --------
    % An integer fold id per observation. fmri_data.predict accepts this as
    % the 'nfolds' argument (a vector => custom holdout sets), and the
    % newapi path reuses the SAME ids via cv_splitter.custom_partition.
    foldvec = mod(0:n-1, nfolds)' + 1;            % deterministic round-robin

    % Synthetic outcomes (same across runs via the seed).
    Ybin  = ones(n, 1); Ybin(2:2:end) = -1;       % balanced +/-1
    Ycont = randn(n, 1);                          % continuous

    rows = {};
    for i = 1:numel(algos)
        spec = algos{i};
        if iscell(spec)
            alg = spec{1};  extra = spec{2};
        else
            alg = spec;     extra = {};
        end
        label = alg;
        if ~isempty(extra), label = [alg ' ' opts_label(extra)]; end

        is_class = strcmp(alg, 'cv_svm');
        if is_class
            dat.Y = Ybin;  task = 'classification';  errtype = 'mcr';
        else
            dat.Y = Ycont; task = 'regression';      errtype = 'mse';
        end

        % --- legacy ---
        [cverr_L, stats_L] = predict(dat, 'algorithm_name', alg, ...
            'nfolds', foldvec, 'error_type', errtype, 'verbose', 0, extra{:});

        % --- newapi (same folds, same error_type, same extra options) ---
        [cverr_N, stats_N] = predict(dat, 'algorithm_name', alg, ...
            'nfolds', foldvec, 'error_type', errtype, 'newapi', 'verbose', 0, extra{:});

        yL = stats_L.yfit(:);  yN = stats_N.yfit(:);
        ok = ~isnan(yL) & ~isnan(yN);
        maxdiff = max(abs(yL(ok) - yN(ok)));
        if sum(ok) > 1 && std(yL(ok)) > 0 && std(yN(ok)) > 0
            ryy = corr(yL(ok), yN(ok));
        else
            ryy = NaN;
        end

        rows(end+1, :) = { label, task, maxdiff, ryy, cverr_L, cverr_N, ...
            abs(cverr_L - cverr_N), maxdiff <= tol }; %#ok<AGROW>

        fprintf(['%-26s [%s]  max|dyfit|=%.3g  corr=%.4f  ' ...
                 'cverr L=%.4f N=%.4f  -> %s\n'], ...
            label, task, maxdiff, ryy, cverr_L, cverr_N, ...
            ternary(maxdiff <= tol, 'IDENTICAL', 'DIFFERS'));
    end

    results = cell2table(rows, 'VariableNames', ...
        {'algorithm','task','max_abs_yfit_diff','corr_yfit', ...
         'cverr_legacy','cverr_newapi','cverr_absdiff','identical'});

    fprintf('\nFolds held constant (%d-fold round-robin). tol=%g.\n', nfolds, tol);
    fprintf(['Expected: cv_svm, cv_pcr, cv_lassopcr (default + lasso_num) IDENTICAL; ' ...
             'cv_svr NEAR-identical (Spider SVR vs fitrsvm solver); ' ...
             'cv_lassopcr estimateparam may differ (new path uses a clean ' ...
             'nested CV — the corrected behavior).\n']);

    % Convenience flags: exact at tol, or near (corr > 0.999).
    results.near_identical = results.corr_yfit > 0.999;
end


function v = ternary(c, a, b)
    if c, v = a; else, v = b; end
end


function s = opts_label(extra)
% Compact label for an extra-options cell, e.g. {'lasso_num',5} -> "(lasso_num=5)".
    parts = {};
    i = 1;
    while i <= numel(extra)
        nm = char(string(extra{i}));
        if i < numel(extra) && ~ischar(extra{i+1}) && ~isstring(extra{i+1})
            parts{end+1} = sprintf('%s=%s', nm, mat2str(extra{i+1})); %#ok<AGROW>
            i = i + 2;
        else
            parts{end+1} = nm; %#ok<AGROW>
            i = i + 1;
        end
    end
    s = ['(' strjoin(parts, ',') ')'];
end
