function pmodel_obj = xval_SVR(X, Y, id, varargin)
%  Thin wrapper: SVR regression using the new @predictive_model API.
%
% :Usage:
% ::
%
%     pmodel_obj = xval_SVR(X, Y, id, varargin)
%
% This is a thin wrapper around `predictive_model` + `crossval`
% (+ optional `grid_search` + `bootstrap`) for SVR regression.
% Subjects are kept together in the same fold via
% `cv_splitter.group_kfold(nfolds)`, and the within-person scoring
% fields (scorediff, crossval_accuracy_within, d_within,
% scores_within_id, Y_within_id, high_vs_low_scores_within_id) are
% populated automatically by `crossval` when `id` has repeated values
% with varying Y.
%
% This implementation REPLACES the legacy ~1100-line xval_SVR with a
% short call to the new pipeline.
%
% :Inputs:
%   X    [n x p] obs x variables predictor matrix
%   Y    [n x 1] continuous outcome vector
%   id   [n x 1] integer grouping vector. Use (1:n)' if no grouping.
%
% :Optional Inputs (name/value or flag form):
%   'modeloptions',     cell    (default {'KernelFunction','linear'})
%   'highdimensional',  logical use fitrlinear (algorithm='linear_svr')
%   'nfolds',           int     (default 10)
%   'nooptimize'        flag    skip hyperparameter search
%   'norepeats'         flag    disable repeated CV
%   'dorepeats',        int     repeats (default 10)
%   'nobootstrap'       flag    skip weight bootstrap
%   'nboot',            int     (default 1000)
%   'noverbose' / 'noplot'      suppress output
%
% :Output:
%   pmodel_obj  @predictive_model. Populated fields include:
%     pm.fitted_values.yfit / .scores / .scorediff /
%                      .scores_within_id / .Y_within_id /
%                      .high_vs_low_scores_within_id
%     pm.weights.w / .boot_w / .z / .p / .fdr_thr / .fdr_sig
%     pm.error_metrics.r2 / .mse / .rmse / .prediction_outcome_r /
%                      .crossval_accuracy_within / .d_within /
%                      .d_singleinterval
%     pm.diagnostics.mult_obs_within_person
%     pm.cv_partition.trIdx / .teIdx / .nfolds
%     pm.ml_model / pm.fold_models / pm.fit_type
%
% :See also:
%   predictive_model, crossval, bootstrap, grid_search, cv_splitter,
%   cv_scorer, xval_SVM, fmri_data.predict
%
% :Example:
% ::
%
%   % synthetic regression on the DPSP Hot-vs-Warm contrast pattern
%   hw_obj = load_image_set('DPSP_hotwarm');
%   X = double(hw_obj.dat');  Y = hw_obj.Y .* randn(size(hw_obj.Y));
%   id = grp2idx(hw_obj.metadata_table.subj_id);
%   pm = xval_SVR(X, Y, id, 'nooptimize', 'norepeats', 'nobootstrap');

    [varargin, dooptimize ] = pop_flag(varargin, 'nooptimize',  true);
    [varargin, dorepeats_d] = pop_flag(varargin, 'norepeats',   1);
    [varargin, dobootstrap] = pop_flag(varargin, 'nobootstrap', true);
    [varargin, verbose    ] = pop_flag(varargin, 'noverbose',   true);
    [varargin, ~          ] = pop_flag(varargin, 'noplot',      true);

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

    if dorepeats_d == 1 && any(strcmpi(p.Parameters, 'dorepeats'))
        dorepeats = p.Results.dorepeats;
    else
        dorepeats = dorepeats_d;
    end

    Y  = Y(:);
    if isempty(id), id = (1:numel(Y))'; end
    id = id(:);

    if p.Results.highdimensional
        algorithm    = 'linear_svr';
        modeloptions = p.Results.modeloptions;
        kf = find(strcmpi(modeloptions, 'KernelFunction'));
        if ~isempty(kf), modeloptions([kf, kf+1]) = []; end
    else
        algorithm    = 'svr';
        modeloptions = p.Results.modeloptions;
    end

    pm = predictive_model( ...
        'algorithm',    algorithm, ...
        'task',         'regression', ...
        'modeloptions', modeloptions, ...
        'random_state', p.Results.random_state);

    % group_kfold (not stratified — Y is continuous)
    splitter = cv_splitter.group_kfold(p.Results.nfolds);
    pm = crossval(pm, X, Y, ...
        'cv',      splitter, ...
        'groups',  id, ...
        'scoring', 'r2');

    if p.Results.dooptimize && strcmp(algorithm, 'svr')
        pg = struct('BoxConstraint', logspace(-3, 3, 5));
        pm = grid_search(pm, X, Y, pg, 'groups', id, 'verbose', verbose);
    end

    if dorepeats > 1
        repeat_scores = nan(dorepeats, 1);
        repeat_scores(1) = pm.error_metrics.r2.value;
        for r = 2:dorepeats
            pm_r = predictive_model( ...
                'algorithm',    pm.algorithm, ...
                'task',         pm.task, ...
                'modeloptions', pm.modeloptions, ...
                'random_state', r);
            pm_r = crossval(pm_r, X, Y, ...
                'cv',      cv_splitter.group_kfold(p.Results.nfolds), ...
                'groups',  id, ...
                'scoring', 'r2');
            repeat_scores(r) = pm_r.error_metrics.r2.value;
        end
        pm.error_metrics.r2_per_repeat = struct( ...
            'value',   repeat_scores, ...
            'descrip', sprintf('%d-repeat CV r2', dorepeats));
    end

    if p.Results.dobootstrap
        pm = bootstrap(pm, X, Y, ...
            'nboot',   p.Results.nboot, ...
            'groups',  id, ...
            'verbose', verbose);
    end

    pmodel_obj = pm;

end


function [vargs, val] = pop_flag(vargs, flag, default)
    val = default;
    hit = find(strcmpi(vargs, flag), 1, 'first');
    if ~isempty(hit)
        vargs(hit) = [];
        if islogical(default)
            val = ~default;
        else
            val = 1;
        end
    end
end
