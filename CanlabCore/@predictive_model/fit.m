function obj = fit(obj, X, Y, varargin)
% fit  Train the model on the full sample (in-sample fit).
%
% :Usage:
% ::
%     pm = predictive_model('algorithm', 'svm', 'task', 'classification');
%     pm = fit(pm, X, Y);
%     pm = fit(pm, X, Y, 'id', subject_id);
%
% Pre-fit pipeline:
%   1. Reseed RNG from obj.random_state for reproducibility.
%   2. Run predictive_model.detect_bad_data(X, Y); remove bad rows /
%      columns; record them in obj.omitted_cases / obj.omitted_features.
%   3. Standardize X if obj.standardize is true; store mean/std in
%      obj.inputParameters so predict() can re-apply the same transform.
%   4. Look up obj.algorithm in predictive_model.algorithm_registry.
%   5. Fit via the registry's fit_fn with defaults + obj.modeloptions.
%   6. Stash trained model in obj.ml_model; extract weights / intercept
%      into obj.weights.{w, intercept}.
%   7. Set obj.fit_type = 'insample'. Append a history line.
%
% :Inputs:
%   X            [n x p] numeric predictor matrix
%   Y            [n x 1] outcome vector
%
% :Optional Inputs (name/value):
%   'id'         [n x 1] grouping vector (e.g. subject id). Stored as
%                obj.id; not used by this fit but available for later
%                crossval.
%
% :Outputs:
%
%   **obj:**
%        the @predictive_model with obj.ml_model (trained MATLAB model),
%        obj.weights.{w, intercept}, obj.Y, obj.omitted_cases /
%        omitted_features, and obj.fit_type = 'insample' populated. Value
%        semantics: capture the return (m = fit(m, X, Y)).
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm');     % 118 imgs (hot/warm), Y=+/-1
%     X = dat.dat'; Y = dat.Y; id = dat.metadata_table.subj_id;
%     pm = predictive_model('algorithm','svm','task','classification');
%     pm = fit(pm, X, Y, 'id', id);             % in-sample fit on all data
%     numel(pm.weights.w)                        % one weight per voxel
%
% :See also:
%   crossval, predict, score, bootstrap

    if isempty(obj.algorithm)
        error('predictive_model:fit:NoAlgorithm', ...
            'Set obj.algorithm (e.g. ''svm'', ''svr'', ''logistic'') before fitting.');
    end

    p = inputParser; p.KeepUnmatched = true;
    addParameter(p, 'id', []);
    parse(p, varargin{:});
    id_in = p.Results.id;

    if ~isempty(obj.random_state), rng(obj.random_state); end

    % 1. Bad-data check (consistent with the xval_* wrappers).
    [oc, of] = predictive_model.detect_bad_data(X, Y);
    if any(oc) || any(of)
        X(oc, :)            = [];
        Y(oc)               = [];
        if ~isempty(id_in), id_in(oc) = []; end
        X(:, of)            = [];
    end
    obj.omitted_cases    = oc;
    obj.omitted_features = of;
    obj.Y                = Y(:);
    if ~isempty(id_in), obj.id = id_in(:); end

    % 2. Standardize if requested.
    if isequal(obj.standardize, true)
        mu = mean(X, 1);
        sd = std(X,  0, 1);
        sd(sd == 0) = 1;
        X = (X - mu) ./ sd;
        if ~isstruct(obj.inputParameters), obj.inputParameters = struct(); end
        obj.inputParameters.standardize_mu = mu;
        obj.inputParameters.standardize_sd = sd;
    end

    % 3a. PCR (principal-components regression) is handled specially: it is
    % not a single MATLAB fit*() object but PCA + OLS on the components,
    % faithful to the legacy fmri_data.predict cv_pcr (and the default
    % cv_lassopcr, which OLS-refits the full model and is identical).
    if strcmpi(obj.algorithm, 'pcr')
        [w, b0] = predictive_model.fit_pcr(X, Y, obj.modeloptions);
        obj.ml_model          = struct('type', 'pcr', 'w', w, 'intercept', b0);
        obj.weights.w         = w;
        obj.weights.intercept = b0;
        if isempty(obj.task), obj.task = 'regression'; end
        obj.fit_type = 'insample';
        obj.history{end+1, 1} = sprintf('fit(pcr): %d obs, %d features', ...
            numel(Y), size(X, 2));
        return
    end

    % 3b. LASSO-PCR: PCA + lasso (shrinkage by 'lasso_num' or nested-CV
    % 'estimateparam') + relaxed-OLS refit on the non-zero components.
    % Faithful to the legacy fmri_data.predict cv_lassopcr.
    if strcmpi(obj.algorithm, 'lassopcr')
        [w, b0, info] = predictive_model.fit_lassopcr(X, Y, obj.modeloptions);
        obj.ml_model          = struct('type', 'pcr', 'w', w, 'intercept', b0);
        obj.weights.w         = w;
        obj.weights.intercept = b0;
        obj.diagnostics.lassopcr = info;
        if isempty(obj.task), obj.task = 'regression'; end
        obj.fit_type = 'insample';
        obj.history{end+1, 1} = sprintf('fit(lassopcr): %d obs, %d features, %d nonzero', ...
            numel(Y), size(X, 2), info.n_nonzero);
        return
    end

    % 3. Look up algorithm.
    reg = predictive_model.algorithm_registry();
    if ~isfield(reg, obj.algorithm)
        error('predictive_model:fit:UnknownAlgorithm', ...
            'Unknown algorithm ''%s''. Registered: %s', ...
            obj.algorithm, strjoin(fieldnames(reg), ', '));
    end
    a = reg.(obj.algorithm);
    if isempty(obj.task)
        obj.task = a.task;
    end

    % 4. Fit.
    opts = [a.defaults, obj.modeloptions];
    mdl = a.fit_fn(X, Y, opts{:});

    % 5. Stash.
    obj.ml_model        = mdl;
    obj.weights.w         = predictive_model.extract_weights(mdl);
    obj.weights.intercept = predictive_model.extract_intercept(mdl);

    obj.fit_type    = 'insample';
    obj.history{end+1, 1} = sprintf('fit(%s): %d obs, %d features', ...
        obj.algorithm, numel(Y), size(X, 2));
end
