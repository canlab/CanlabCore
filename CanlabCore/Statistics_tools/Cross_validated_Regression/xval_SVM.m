function pmodel_obj = xval_SVM(X, Y, id, varargin)
%  Thin wrapper: SVM classification using the new @predictive_model API.
%
% :Usage:
% ::
%
%     pmodel_obj = xval_SVM(X, Y, id, varargin)
%
% This is a thin wrapper around `predictive_model` + `crossval`
% (+ optional `grid_search` + `bootstrap`) for binary SVM
% classification. Subjects are kept together in the same fold via
% `cv_splitter.stratified_group_kfold(nfolds)`, and the within-person
% scoring fields (`scorediff`, `crossval_accuracy_within`, `d_within`,
% `scores_within_id`, `Y_within_id`) are populated automatically by
% `crossval` when `id` has repeated values with varying Y.
%
% This implementation REPLACES the legacy ~1350-line xval_SVM with a
% short call to the new pipeline. The output predictive_model carries
% the same fields users have been reading (yfit, w, crossval_accuracy,
% d_singleinterval, d_within, dist_from_hyperplane_xval, fdr_sig, etc.)
% so existing downstream code keeps working.
%
% :Inputs:
%   X    [n x p] obs x variables predictor matrix
%   Y    [n x 1] outcome vector, effects-coded (+1, -1)
%   id   [n x 1] integer grouping vector (e.g. participant id). Use
%        (1:n)' if there is no grouping.
%
% :Optional Inputs (name/value or flag form):
%   'modeloptions',     cell        forwarded to fitcsvm/fitclinear
%                                   (default {'KernelFunction','linear'})
%   'highdimensional',  logical     if true, use fitclinear (algorithm =
%                                   'linear_svm'); recommended for wide
%                                   data (p > a few thousand)
%   'nfolds',           positive int (default 10)
%   'nooptimize'        flag         skip hyperparameter grid_search
%                                    (default: do optimize via a small
%                                    BoxConstraint grid; legacy xval_SVM
%                                    used Bayesian fitcsvm
%                                    OptimizeHyperparameters — the new
%                                    wrapper does grid_search instead)
%   'norepeats'         flag         disable repeated cross-validation
%   'dorepeats',        int          number of CV repeats (default 10)
%   'nobootstrap'       flag         skip bootstrap of weights
%   'nboot',            int          (default 1000)
%   'noverbose' / 'noplot'           suppress output (passed through)
%
% :Output:
%   pmodel_obj   a @predictive_model object. Same field layout as the
%                legacy xval_SVM:
%                   pm.fitted_values.yfit / .scores / .scorediff /
%                                    .dist_from_hyperplane_xval /
%                                    .scores_within_id / .Y_within_id
%                   pm.weights.w / .boot_w / .z / .p / .fdr_thr /
%                                  .fdr_sig / .thresh_fdr
%                   pm.error_metrics.crossval_accuracy / .d_singleinterval /
%                                    .crossval_accuracy_within / .d_within /
%                                    .balanced_accuracy / ...
%                   pm.diagnostics.mult_obs_within_person
%                   pm.cv_partition.trIdx / .teIdx / .nfolds
%                   pm.ml_model        the full-sample trained model
%                   pm.fold_models     per-fold trained models
%                   pm.fit_type        'crossval'
%                   pm.history         provenance trail
%
% :See also:
%   predictive_model, weight_map_object, crossval, bootstrap, grid_search,
%   cv_splitter, cv_scorer, fmri_data.fit_predictive_model, fmri_data.predict
%
% :Example:
% ::
%
%   hw_obj = load_image_set('DPSP_hotwarm');
%   X  = double(hw_obj.dat');
%   Y  = hw_obj.Y;
%   id = grp2idx(hw_obj.metadata_table.subj_id);
%   pm = xval_SVM(X, Y, id, 'highdimensional', true, ...
%                'nooptimize', 'norepeats', 'nobootstrap');
%   pm.error_metrics.crossval_accuracy.value          % cv accuracy
%   pm.error_metrics.crossval_accuracy_within.value   % within-person acc
%   pm.error_metrics.d_within.value                   % effect size
%   % This wrapper trains on matrices and never sees an image object, so it
%   % cannot build a weight map. Attach one with any reference image in the
%   % same space, then montage(pm) / surface(pm) work with no source:
%   pm = weight_map_object(pm, hw_obj);               % cache weight @statistic_image
%   montage(pm); surface(pm);
%   [~, si] = weight_map_object(pm, hw_obj);          % (or get the image directly)

    % ----- parse flag-style options that aren't name/value -----
    [varargin, dooptimize ] = pop_flag(varargin, 'nooptimize',  true);
    [varargin, dorepeats_d] = pop_flag(varargin, 'norepeats',   1);    % 1 = no repeat
    [varargin, dobootstrap] = pop_flag(varargin, 'nobootstrap', true);
    [varargin, verbose    ] = pop_flag(varargin, 'noverbose',   true);
    [varargin, ~          ] = pop_flag(varargin, 'noplot',      true); % accepted but unused

    p = inputParser; p.KeepUnmatched = true;
    addParameter(p, 'modeloptions',    {'KernelFunction', 'linear'});
    addParameter(p, 'highdimensional', false);
    addParameter(p, 'nfolds',          10);
    addParameter(p, 'dooptimize',      dooptimize);
    addParameter(p, 'dorepeats',       dorepeats_d);
    addParameter(p, 'dobootstrap',     dobootstrap);
    addParameter(p, 'nboot',           1000);
    addParameter(p, 'random_state',    []);
    parse(p, varargin{:});

    % If the user passed 'dorepeats' explicitly, that takes precedence
    % over the 'norepeats' flag default (1).
    if dorepeats_d == 1 && any(strcmpi(p.Parameters, 'dorepeats'))
        dorepeats = p.Results.dorepeats;
    else
        dorepeats = dorepeats_d;
    end
    % Default dorepeats when user passed nothing: 1 (single CV pass).

    Y  = Y(:);
    if isempty(id), id = (1:numel(Y))'; end
    id = id(:);

    % ----- build the predictive_model -----
    if p.Results.highdimensional
        algorithm    = 'linear_svm';
        modeloptions = p.Results.modeloptions;
        % fitclinear doesn't take KernelFunction; strip if present
        kf = find(strcmpi(modeloptions, 'KernelFunction'));
        if ~isempty(kf), modeloptions([kf, kf+1]) = []; end
    else
        algorithm    = 'svm';
        modeloptions = p.Results.modeloptions;
    end

    pm = predictive_model( ...
        'algorithm',    algorithm, ...
        'task',         'classification', ...
        'modeloptions', modeloptions, ...
        'random_state', p.Results.random_state);

    % ----- cross-validation (with stratified-group-kfold for the
    % canonical "no subject spans train/test" guarantee) -----
    splitter = cv_splitter.stratified_group_kfold(p.Results.nfolds);
    pm = crossval(pm, X, Y, ...
        'cv',      splitter, ...
        'groups',  id, ...
        'scoring', 'balanced_accuracy');

    % ----- optional grid_search over BoxConstraint (kernel SVM only;
    % fitclinear's tuning is different so we skip then) -----
    if p.Results.dooptimize && strcmp(algorithm, 'svm')
        pg = struct('BoxConstraint', logspace(-3, 3, 5));
        pm = grid_search(pm, X, Y, pg, 'groups', id, 'verbose', verbose);
    end

    % ----- optional repeated CV (aggregate per-repeat scores) -----
    if dorepeats > 1
        repeat_scores = nan(dorepeats, 1);
        repeat_scores(1) = pm.error_metrics.balanced_accuracy.value;
        for r = 2:dorepeats
            pm_r = predictive_model( ...
                'algorithm',    pm.algorithm, ...
                'task',         pm.task, ...
                'modeloptions', pm.modeloptions, ...
                'random_state', r);          % different seed per repeat
            pm_r = crossval(pm_r, X, Y, ...
                'cv',      cv_splitter.stratified_group_kfold(p.Results.nfolds), ...
                'groups',  id, ...
                'scoring', 'balanced_accuracy');
            repeat_scores(r) = pm_r.error_metrics.balanced_accuracy.value;
        end
        pm.error_metrics.balanced_accuracy_per_repeat = struct( ...
            'value',   repeat_scores, ...
            'descrip', sprintf('%d-repeat CV balanced accuracies', dorepeats));
    end

    % ----- optional bootstrap -----
    if p.Results.dobootstrap
        pm = bootstrap(pm, X, Y, ...
            'nboot',   p.Results.nboot, ...
            'groups',  id, ...
            'verbose', verbose);
    end

    % (crossval auto-populates pm.error_metrics.crossval_accuracy as a
    % 100x scaled alias of balanced_accuracy for classification, so we
    % don't need to add it here.)

    pmodel_obj = pm;

end


% ----- local helper: detect & remove a flag-style option in varargin
function [vargs, val] = pop_flag(vargs, flag, default)
    val = default;
    hit = find(strcmpi(vargs, flag), 1, 'first');
    if ~isempty(hit)
        vargs(hit) = [];
        % flags invert their default: 'nooptimize' turns off
        % optimization, etc.
        if islogical(default)
            val = ~default;
        else
            val = 1;   % numeric default (e.g. dorepeats=1 means no repeats)
        end
    end
end
