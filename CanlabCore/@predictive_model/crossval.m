function obj = crossval(obj, X, Y, varargin)
% crossval  Cross-validate the model: refit per fold, predict held-out, score.
%
% :Usage:
% ::
%     pm = predictive_model('algorithm','svm','task','classification');
%     pm = crossval(pm, X, Y, 'cv', cv_splitter.stratified_kfold(5));
%     pm = crossval(pm, X, Y, 'groups', subject_id, ...
%                              'cv',      cv_splitter.group_kfold(10), ...
%                              'scoring', 'roc_auc');
%
% Pipeline:
%   1. Pre-fit bad-data check (same as fit()); drop rows/columns.
%   2. Pick cv splitter (obj.cv) and scorer (obj.scorer); default to
%      stratified_kfold(5) + balanced_accuracy for classification,
%      kfold(5) + r2 for regression.
%   3. For each fold:
%        - clone(obj); fit on training rows; predict held-out;
%        - aggregate yfit and continuous scores;
%        - record per-fold score and fit time.
%   4. Re-fit on full data so obj.ml_model / obj.weights are the
%      canonical "ship-it" model.
%   5. Populate:
%        obj.fitted_values.yfit / .scores / .score_type
%        obj.cv_partition.{trIdx, teIdx, nfolds}
%        obj.fold_models     (per-fold ml_models, or full pm clones
%                             if 'return_estimator' is true)
%        obj.error_metrics.<scorer_name>            mean across folds
%        obj.error_metrics.<scorer_name>_per_fold   per-fold vector
%        obj.error_metrics.fit_time_per_fold        timings
%        obj.fit_type = 'crossval'
%
% :Inputs:
%
%   **obj:**
%        a @predictive_model with obj.algorithm set.
%
%   **X:**
%        [n x p] predictor matrix.
%
%   **Y:**
%        [n x 1] outcome vector.
%
% :Optional Inputs (name/value):
%   'groups'           grouping vector (e.g. subject id) for grouped splitters
%   'cv'               cv_splitter object; if omitted, uses obj.cv or default
%   'scoring'          scorer name (e.g. 'roc_auc') OR a cv_scorer object;
%                      if omitted, uses obj.scorer or default-for-task
%   'return_estimator' (default false) — store full pm clones in fold_models
%                      instead of just the inner ml_models
%
% :Outputs:
%
%   **obj:**
%        the @predictive_model with cross-validated state populated:
%        obj.fitted_values.{yfit, scores}, obj.error_metrics.<scorer>
%        (.value/.descrip tuples) plus per-fold vectors, obj.cv_partition,
%        obj.fold_models, and obj.fit_type = 'crossval'. obj.ml_model /
%        obj.weights are refit on the FULL data (the ship-it model). When
%        'groups' is given, within-person forced-choice stats are also
%        added (see compute_within_person_stats).
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm');
%     X = dat.dat'; Y = dat.Y; id = dat.metadata_table.subj_id;
%     pm = predictive_model('algorithm','svm','task','classification');
%     % leave-whole-subject-out, stratified on Y within the held-out folds
%     pm = crossval(pm, X, Y, 'groups', id, ...
%                   'cv', cv_splitter.stratified_group_kfold(5));
%     pm.error_metrics.crossval_accuracy.value      % cv accuracy (%)
%     pm.error_metrics.crossval_accuracy_within.value  % forced-choice acc
%
% :See also:
%   fit, predict, score, bootstrap, cv_splitter, cv_scorer,
%   report_accuracy, summary

    p = inputParser; p.KeepUnmatched = true;
    addParameter(p, 'groups',           []);
    addParameter(p, 'cv',               []);
    addParameter(p, 'scoring',          []);
    addParameter(p, 'return_estimator', false);
    parse(p, varargin{:});

    groups            = p.Results.groups;
    return_estimator  = p.Results.return_estimator;

    % Normalize non-numeric groups (cell of subject-id strings, categorical,
    % string array) to integer labels so within-person stats and the
    % splitters can compare them.
    if ~isempty(groups) && ~(isnumeric(groups) || islogical(groups))
        [~, ~, groups] = unique(groups(:), 'stable');
    end

    if ~isempty(p.Results.cv), obj.cv = p.Results.cv; end
    if ~isempty(p.Results.scoring)
        if isa(p.Results.scoring, 'cv_scorer')
            obj.scorer = p.Results.scoring;
        else
            obj.scorer = cv_scorer.make(p.Results.scoring);
        end
    end

    if isempty(obj.algorithm)
        error('predictive_model:crossval:NoAlgorithm', ...
            'Set obj.algorithm before crossval.');
    end

    % Pre-fit bad-data check (same as fit()).
    [oc, of] = predictive_model.detect_bad_data(X, Y);
    if any(oc) || any(of)
        X(oc, :) = []; Y(oc) = []; X(:, of) = [];
        if ~isempty(groups), groups(oc) = []; end
    end

    Y = Y(:);
    n = numel(Y);

    % Infer task if not set.
    if isempty(obj.task)
        if numel(unique(Y(~isnan(Y)))) <= 2
            obj.task = 'classification';
        else
            obj.task = 'regression';
        end
    end

    % Pick cv splitter.
    if isempty(obj.cv)
        if strcmp(obj.task, 'classification')
            obj.cv = cv_splitter.stratified_kfold(5);
        else
            obj.cv = cv_splitter.kfold(5);
        end
    end
    if ~isempty(obj.random_state) && isempty(obj.cv.random_state)
        obj.cv.random_state = obj.random_state;
    end

    % Pick scorer.
    if isempty(obj.scorer)
        obj.scorer = cv_scorer.default_for_task(obj.task);
    end

    splits = obj.cv.split(X, Y, groups);
    nfolds = numel(splits);

    yfit_cv         = nan(n, 1);
    scores_cv       = nan(n, 1);
    per_fold_score  = nan(nfolds, 1);
    fit_time        = zeros(nfolds, 1);
    fold_estimators = cell(1, nfolds);

    needs_cont = obj.scorer.needs_continuous;

    for f = 1:nfolds
        tr = splits(f).trIdx;
        te = splits(f).teIdx;
        if ~any(tr) || ~any(te)
            continue
        end

        m = clone(obj);
        t0 = tic;
        m = fit(m, X(tr, :), Y(tr));
        fit_time(f) = toc(t0);

        if needs_cont
            [yf, sf] = predict(m, X(te, :));
            yfit_cv(te)   = yf;
            % Use the rightmost score column (class-1 / regression yhat).
            if size(sf, 2) >= 1
                scores_cv(te) = sf(:, end);
            end
            per_fold_score(f) = obj.scorer.score(Y(te), yf, sf);
        else
            [yf, sf] = predict(m, X(te, :));
            yfit_cv(te) = yf;
            if size(sf, 2) >= 1
                scores_cv(te) = sf(:, end);
            end
            per_fold_score(f) = obj.scorer.score(Y(te), yf);
        end

        if return_estimator
            fold_estimators{f} = m;
        else
            fold_estimators{f} = m.ml_model;
        end
    end

    % Final: refit on full data so the canonical obj.ml_model / .weights.w
    % is the "ship-it" model, not a per-fold model.
    obj = fit(obj, X, Y, 'id', groups);

    % Store cross-validation outputs (overwriting the in-sample fit_type
    % set by fit()).
    obj.fitted_values.yfit       = yfit_cv;
    if any(~isnan(scores_cv))
        obj.fitted_values.scores     = scores_cv;
        obj.fitted_values.score_type = ternary(needs_cont, 'probability', 'distance');
        % Legacy alias for xval_SVM-style consumers: dist_from_hyperplane_xval
        % is the same data, named for the SVM use case.
        if strcmpi(obj.task, 'classification')
            obj.fitted_values.dist_from_hyperplane_xval = scores_cv;
        end
    end
    obj.fit_type = 'crossval';

    % Within-person scoring (auto when groups is provided and any group
    % has >=2 obs with varying Y). Populates fitted_values.scorediff /
    % scores_within_id / Y_within_id and error_metrics.{
    % crossval_accuracy_within, d_within, d_singleinterval}, plus a
    % diagnostics.mult_obs_within_person flag.
    if ~isempty(groups)
        wps = predictive_model.compute_within_person_stats( ...
            Y, scores_cv, groups, obj.task);
        obj.fitted_values.scorediff             = wps.scorediff;
        obj.fitted_values.scores_within_id      = wps.scores_within_id;
        obj.fitted_values.Y_within_id           = wps.Y_within_id;
        if ~isempty(wps.high_vs_low_scores_within_id)
            obj.fitted_values.high_vs_low_scores_within_id = ...
                wps.high_vs_low_scores_within_id;
        end
        obj.diagnostics.mult_obs_within_person  = wps.mult_obs_within_person;
        if ~isnan(wps.crossval_accuracy_within)
            obj.error_metrics.crossval_accuracy_within = struct( ...
                'value',   wps.crossval_accuracy_within, ...
                'descrip', 'forced-choice accuracy: fraction of subjects with scorediff > 0');
        end
        if ~isnan(wps.d_within)
            obj.error_metrics.d_within = struct( ...
                'value',   wps.d_within, ...
                'descrip', 'Cohen''s d on paired within-subject scorediff');
        end
        if ~isnan(wps.d_singleinterval)
            obj.error_metrics.d_singleinterval = struct( ...
                'value',   wps.d_singleinterval, ...
                'descrip', 'pooled standardized mean diff of cv scores between classes');
        end
    end

    obj.fold_models               = fold_estimators;
    obj.cv_partition.nfolds       = nfolds;
    obj.cv_partition.trIdx        = arrayfun(@(s) s.trIdx, splits, 'uniform', false);
    obj.cv_partition.teIdx        = arrayfun(@(s) s.teIdx, splits, 'uniform', false);

    obj.error_metrics.(obj.scorer.name) = struct( ...
        'value',   mean(per_fold_score, 'omitnan'), ...
        'descrip', sprintf('mean %s across %d folds (greater_is_better=%d)', ...
                           obj.scorer.name, nfolds, obj.scorer.greater_is_better));

    % Legacy alias: xval_SVM/SVR consumers read pm.error_metrics.crossval_accuracy.
    % For classification scorers in [0,1], expose the 0-100 scaled value
    % under that name. (.value as percent matches legacy CANlab convention.)
    if strcmpi(obj.task, 'classification') ...
            && ~isfield(obj.error_metrics, 'crossval_accuracy') ...
            && any(strcmpi(obj.scorer.name, {'accuracy', 'balanced_accuracy'}))
        obj.error_metrics.crossval_accuracy = struct( ...
            'value',   100 * obj.error_metrics.(obj.scorer.name).value, ...
            'descrip', sprintf('cv accuracy as a percent (= 100 * %s)', obj.scorer.name));
    end

    % Regression out-of-sample performance from the pooled cross-validated
    % predictions: Pearson r (xval_SVR consumers read prediction_outcome_r)
    % plus the two section-39.4 R^2 variants (predicted_r2, out_of_sample_r2).
    % These pool all held-out predictions rather than averaging per-fold r2.
    if strcmpi(obj.task, 'regression') ...
            && ~isempty(yfit_cv) && any(~isnan(yfit_cv))
        rm = predictive_model.regression_metrics_from_cv( ...
            Y, yfit_cv, obj.cv_partition.trIdx, obj.cv_partition.teIdx);

        if ~isfield(obj.error_metrics, 'prediction_outcome_r') && ~isnan(rm.prediction_outcome_r)
            obj.error_metrics.prediction_outcome_r = struct( ...
                'value',   rm.prediction_outcome_r, ...
                'descrip', 'Pearson r between cv predictions and Y');
        end
        if ~isnan(rm.predicted_r2)
            obj.error_metrics.predicted_r2 = struct( ...
                'value',   rm.predicted_r2, ...
                'descrip', 'predicted R^2 = 1 - PRESS/SST (grand-mean denominator; Wager & Lindquist Ch.39.4)');
        end
        if ~isnan(rm.out_of_sample_r2)
            obj.error_metrics.out_of_sample_r2 = struct( ...
                'value',   rm.out_of_sample_r2, ...
                'descrip', 'out-of-sample R^2 = 1 - PRESS / SS(y - per-fold train mean) (Ch.39.4)');
        end
    end
    obj.error_metrics.([obj.scorer.name '_per_fold']) = struct( ...
        'value',   per_fold_score, ...
        'descrip', sprintf('per-fold %s', obj.scorer.name));
    obj.error_metrics.fit_time_per_fold = struct( ...
        'value', fit_time, 'descrip', 'seconds per fold');

    obj.history{end+1, 1} = sprintf('crossval(%s, %s, scorer=%s): mean=%.3f', ...
        obj.algorithm, obj.cv.type, obj.scorer.name, ...
        mean(per_fold_score, 'omitnan'));
end


% -------- small helper (avoids if/else nesting in struct() args) ----------
function v = ternary(cond, a, b)
    if cond, v = a; else, v = b; end
end
