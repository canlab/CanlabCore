function obj = crossval(obj, X, Y, varargin)
% crossval  Leakage-free cross-validation of a pipeline.
%
% Every transform step is refit on the TRAINING rows of each fold (the
% held-out rows never touch the transformer fit), then the estimator is
% fit on the transformed training data and applied to the transformed
% held-out data. This is the correct way to cross-validate
% PCA-then-regression, standardize-then-SVM, select-then-classify, etc.
%
% :Usage:
% ::
%     est  = predictive_model('algorithm','ridge','task','regression');
%     pipe = pipeline({ {'pca','k',20} }, est);
%     pipe = crossval(pipe, X, Y);
%     pipe = crossval(pipe, X, Y, 'cv', cv_splitter.kfold(10), 'groups', id);
%
% :Inputs:
%
%   **obj:**  a @pipeline (with obj.estimator set).
%   **X:**    [n x p] predictor matrix.
%   **Y:**    [n x 1] outcome vector.
%
% :Optional Inputs (name/value):
%
%   **'cv':**       a cv_splitter; defaults to obj.cv, else stratified_kfold(5)
%                   for classification / kfold(5) for regression.
%   **'groups':**   grouping vector for grouped splitters (e.g. subject id).
%   **'scoring':**  a cv_scorer or its name; defaults to the task default.
%
% :Outputs:
%
%   Populates obj.fitted_values.yfit / .scores, obj.error_metrics
%   (<scorer> mean + per-fold), obj.cv_partition, obj.fit_type='crossval',
%   and refits all steps + estimator on the full data so obj is the
%   ship-it model. obj.weights.w holds the back-projected input-space
%   weights of that full-data fit.
%
% :See also:
%   pipeline, predictive_model.crossval, cv_splitter, cv_scorer

    if isempty(obj.estimator)
        error('pipeline:crossval:NoEstimator', 'Set obj.estimator before crossval.');
    end

    p = inputParser; p.KeepUnmatched = true;
    addParameter(p, 'groups',  []);
    addParameter(p, 'cv',      []);
    addParameter(p, 'scoring', []);
    parse(p, varargin{:});
    groups = p.Results.groups;

    if ~isempty(p.Results.cv), obj.cv = p.Results.cv; end
    if ~isempty(p.Results.scoring)
        if isa(p.Results.scoring, 'cv_scorer')
            obj.scorer = p.Results.scoring;
        else
            obj.scorer = cv_scorer.make(p.Results.scoring);
        end
    end

    % Infer task from estimator (fall back to Y).
    task = obj.estimator.task;
    if isempty(task)
        if numel(unique(Y(~isnan(Y)))) <= 2, task = 'classification'; else, task = 'regression'; end
        obj.estimator.task = task;
    end

    % Bad-CASE check only (drop NaN/Inf rows). We deliberately do NOT drop
    % input feature columns here: the transform steps tolerate constant /
    % zero-variance columns (zscore floors SD, PCA via SVD is rank-safe),
    % and keeping the full input width means back-projected weights stay
    % aligned to the source voxel space.
    oc = predictive_model.detect_bad_data(X, Y);
    if any(oc)
        X(oc, :) = []; Y(oc) = [];
        if ~isempty(groups), groups(oc) = []; end
    end
    Y = Y(:); n = numel(Y);

    % Splitter.
    if isempty(obj.cv)
        if strcmp(task, 'classification')
            obj.cv = cv_splitter.stratified_kfold(5);
        else
            obj.cv = cv_splitter.kfold(5);
        end
    end
    if ~isempty(obj.random_state) && isempty(obj.cv.random_state)
        obj.cv.random_state = obj.random_state;
    end

    % Scorer.
    if isempty(obj.scorer), obj.scorer = cv_scorer.default_for_task(task); end
    needs_cont = obj.scorer.needs_continuous;

    splits = obj.cv.split(X, Y, groups);
    nfolds = numel(splits);

    yfit_cv        = nan(n, 1);
    scores_cv      = nan(n, 1);
    per_fold_score = nan(nfolds, 1);

    for f = 1:nfolds
        tr = splits(f).trIdx; te = splits(f).teIdx;
        if ~any(tr) || ~any(te), continue; end

        m = clone(obj);                 % fresh steps + cloned estimator
        m = fit(m, X(tr, :), Y(tr));    % refits EVERY step on train only
        [yf, sf] = predict(m, X(te, :));
        yfit_cv(te) = yf;
        if ~isempty(sf), scores_cv(te) = sf(:, end); end
        if needs_cont
            per_fold_score(f) = obj.scorer.score(Y(te), yf, sf);
        else
            per_fold_score(f) = obj.scorer.score(Y(te), yf);
        end
    end

    % Refit on the full data (ship-it model + back-projected weights).
    obj = fit(obj, X, Y);

    obj.fitted_values.yfit = yfit_cv;
    if any(~isnan(scores_cv))
        obj.fitted_values.scores     = scores_cv;
        obj.fitted_values.score_type = ternary(needs_cont, 'probability', 'distance');
        if strcmp(task, 'classification')
            obj.fitted_values.dist_from_hyperplane_xval = scores_cv;
        end
    end
    obj.fit_type = 'crossval';

    obj.cv_partition.nfolds = nfolds;
    obj.cv_partition.trIdx  = arrayfun(@(s) s.trIdx, splits, 'uniform', false);
    obj.cv_partition.teIdx  = arrayfun(@(s) s.teIdx, splits, 'uniform', false);

    obj.error_metrics.(obj.scorer.name) = struct( ...
        'value',   mean(per_fold_score, 'omitnan'), ...
        'descrip', sprintf('mean %s across %d folds', obj.scorer.name, nfolds));
    obj.error_metrics.([obj.scorer.name '_per_fold']) = struct( ...
        'value', per_fold_score, 'descrip', sprintf('per-fold %s', obj.scorer.name));

    step_names = cellfun(@(s) s.name, obj.steps, 'uniform', false);
    obj.history{end+1, 1} = sprintf('crossval [%s] -> %s (%s): mean %s = %.3f', ...
        strjoin(step_names, ' | '), char(string(obj.estimator.algorithm)), ...
        obj.cv.type, obj.scorer.name, mean(per_fold_score, 'omitnan'));
end


function v = ternary(c, a, b)
    if c, v = a; else, v = b; end
end
