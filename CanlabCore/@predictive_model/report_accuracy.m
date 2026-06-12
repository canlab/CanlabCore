function acc = report_accuracy(obj, varargin)
% report_accuracy  Model-type-relevant performance metrics, as a struct + table.
%
% Collects the performance metrics that matter for the model's task into a
% single struct and (by default) prints them as a labelled table. For
% classification it pulls discrimination/threshold metrics from the
% cross-validated ROC (AUC, sensitivity, specificity, PPV, NPV) plus
% accuracy / balanced accuracy and forced-choice accuracy when available;
% for regression it reports the prediction-outcome correlation and the
% out-of-sample R-squared variants from Wager & Lindquist Ch.39.4
% (predicted_r2 = 1 - PRESS/SST, and out_of_sample_r2), alongside RMSE / MAE
% / MSE.
%
% Values are read from obj.error_metrics when present and otherwise computed
% from obj.fitted_values (cross-validated predictions / scores) and obj.Y, so
% the report reflects out-of-sample performance whenever the model was
% cross-validated.
%
% :Usage:
% ::
%     report_accuracy(pm);             % print the table
%     acc = report_accuracy(pm);       % struct of metrics (still prints)
%     acc = report_accuracy(pm, 'noverbose');   % struct only, no printing
%
% :Inputs:
%
%   **obj:**
%        a fitted/cross-validated @predictive_model.
%
% :Optional Inputs (name/value):
%
%   **'doverbose' / 'verbose':**
%        logical, default true. Print the formatted table. 'noverbose' is a
%        shorthand for 'doverbose', false.
%
% :Outputs:
%
%   **acc:**
%        struct of numeric metrics (NaN where not computable). Always has a
%        .task field; the remaining fields depend on the task:
%          classification: accuracy, balanced_accuracy, auc, sensitivity,
%                          specificity, ppv, npv, d, forced_choice_accuracy,
%                          accuracy_p, n
%          regression:     prediction_outcome_r, predicted_r2,
%                          out_of_sample_r2, rmse, mae, mse, n
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm', 'noverbose');
%     X = dat.dat'; Y = dat.Y; id = dat.metadata_table.subj_id;
%     pm  = predictive_model('algorithm','svm','task','classification');
%     pm  = crossval(pm, X, Y, 'groups', id);
%     report_accuracy(pm);
%
%     bmrk3 = load_image_set('bmrk3', 'noverbose');
%     [~,~,sid] = unique(bmrk3.metadata_table.subject_id, 'stable');
%     [~,~,~,pm_r] = predict(bmrk3, 'algorithm_name','cv_pcr', ...
%                            'nfolds', mod(sid,5)+1, 'newapi');
%     acc = report_accuracy(pm_r);
%     fprintf('predicted R^2 = %.3f\n', acc.predicted_r2);
%
% :See also:
%   summary, rocplot, crossval, predictive_model.regression_metrics_from_cv

    p = inputParser; p.KeepUnmatched = true;
    addParameter(p, 'doverbose', true, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'verbose',   [],   @(x) islogical(x) || isnumeric(x));
    % 'noverbose' flag form
    keep = true(size(varargin)); noverbose = false;
    for i = 1:numel(varargin)
        if (ischar(varargin{i}) || isstring(varargin{i})) && strcmpi(varargin{i}, 'noverbose')
            noverbose = true; keep(i) = false;
        end
    end
    parse(p, varargin{keep});
    doverbose = logical(p.Results.doverbose);
    if ~isempty(p.Results.verbose), doverbose = logical(p.Results.verbose); end
    if noverbose, doverbose = false; end

    em = obj.error_metrics;
    is_clf = is_classifier(obj);

    if is_clf
        acc = local_classification(obj, em);
    else
        acc = local_regression(obj, em);
    end

    if doverbose
        local_print(obj, acc);
    end
end


% -------------------------------------------------------------------------
function acc = local_classification(obj, em)
    acc = struct('task', 'classification', ...
        'accuracy', NaN, 'balanced_accuracy', NaN, 'auc', NaN, ...
        'sensitivity', NaN, 'specificity', NaN, 'ppv', NaN, 'npv', NaN, ...
        'd', NaN, 'forced_choice_accuracy', NaN, 'accuracy_p', NaN, 'n', NaN);

    acc.balanced_accuracy = first_num(predictive_model.metric_value(em, 'balanced_accuracy'));
    acc.forced_choice_accuracy = first_num(predictive_model.metric_value(em, 'crossval_accuracy_within'));
    acc.d = first_num(predictive_model.metric_value(em, 'd_within'));
    if isnan(acc.d), acc.d = first_num(predictive_model.metric_value(em, 'd_singleinterval')); end

    % Threshold + discrimination metrics from the cross-validated ROC.
    % Capture roc_plot's own console output so the table stays clean.
    try
        [~] = evalc('ROC = rocplot(obj, ''noplot'');');
        acc.auc         = first_num(ROC.AUC);
        acc.sensitivity = first_num(ROC.sensitivity);
        acc.specificity = first_num(ROC.specificity);
        acc.ppv         = first_num(ROC.PPV);
        acc.accuracy    = first_num(ROC.accuracy);
        if isfield(ROC, 'accuracy_p'), acc.accuracy_p = first_num(ROC.accuracy_p); end
        if isfield(ROC, 'N'),          acc.n          = first_num(ROC.N); end
        if isfield(ROC, 'observations')
            tn = sum(ROC.observations.trueneg);
            fn = sum(ROC.observations.falseneg);
            if (tn + fn) > 0, acc.npv = tn / (tn + fn); end
        end
    catch
        % no continuous scores / not binary — leave ROC-derived fields NaN
    end

    if isnan(acc.accuracy)
        acc.accuracy = first_num(predictive_model.metric_value(em, 'accuracy'));
    end
    if isnan(acc.n) && ~isempty(obj.Y), acc.n = sum(~isnan(obj.Y)); end
end


% -------------------------------------------------------------------------
function acc = local_regression(obj, em)
    acc = struct('task', 'regression', ...
        'prediction_outcome_r', NaN, 'predicted_r2', NaN, ...
        'out_of_sample_r2', NaN, 'rmse', NaN, 'mae', NaN, 'mse', NaN, 'n', NaN);

    acc.prediction_outcome_r = first_num(predictive_model.metric_value(em, 'prediction_outcome_r'));
    acc.predicted_r2     = first_num(predictive_model.metric_value(em, 'predicted_r2'));
    acc.out_of_sample_r2 = first_num(predictive_model.metric_value(em, 'out_of_sample_r2'));
    acc.rmse = first_num(predictive_model.metric_value(em, 'rmse'));
    acc.mae  = first_num(predictive_model.metric_value(em, 'mae'));
    acc.mse  = first_num(predictive_model.metric_value(em, 'mse'));

    % Fill any missing from the pooled cross-validated predictions.
    yfit = [];
    if isstruct(obj.fitted_values) && isfield(obj.fitted_values, 'yfit')
        yfit = obj.fitted_values.yfit;
    end
    need = any(isnan([acc.prediction_outcome_r acc.predicted_r2 ...
                      acc.out_of_sample_r2 acc.rmse acc.mae acc.mse]));
    if need && ~isempty(yfit) && ~isempty(obj.Y)
        tr = {}; te = {};
        if isstruct(obj.cv_partition)
            if isfield(obj.cv_partition, 'trIdx'), tr = obj.cv_partition.trIdx; end
            if isfield(obj.cv_partition, 'teIdx'), te = obj.cv_partition.teIdx; end
        end
        rm = predictive_model.regression_metrics_from_cv(obj.Y, yfit, tr, te);
        acc.prediction_outcome_r = fill(acc.prediction_outcome_r, rm.prediction_outcome_r);
        acc.predicted_r2     = fill(acc.predicted_r2,     rm.predicted_r2);
        acc.out_of_sample_r2 = fill(acc.out_of_sample_r2, rm.out_of_sample_r2);
        acc.rmse = fill(acc.rmse, rm.rmse);
        acc.mae  = fill(acc.mae,  rm.mae);
        acc.mse  = fill(acc.mse,  rm.mse);
    end
    if ~isempty(obj.Y), acc.n = sum(~isnan(obj.Y)); end
end


% -------------------------------------------------------------------------
function local_print(obj, acc)
    src = ternary(strcmpi(char(string(obj.fit_type)), 'crossval'), ...
                  'cross-validated', char(string(obj.fit_type)));
    fprintf('\n  Performance (%s, task=%s)\n', src, acc.task);
    fprintf('  %s\n', repmat('-', 1, 44));
    if strcmpi(acc.task, 'classification')
        rows = { ...
            'Accuracy',               acc.accuracy,               '%.1f%%',  100; ...
            'Balanced accuracy',      acc.balanced_accuracy,      '%.1f%%',  100; ...
            'Forced-choice accuracy', acc.forced_choice_accuracy, '%.1f%%',  1;   ... % already percent
            'AUC',                    acc.auc,                    '%.3f',    1;   ...
            'Sensitivity',            acc.sensitivity,            '%.3f',    1;   ...
            'Specificity',            acc.specificity,            '%.3f',    1;   ...
            'PPV',                    acc.ppv,                    '%.3f',    1;   ...
            'NPV',                    acc.npv,                    '%.3f',    1;   ...
            'Effect size (d)',        acc.d,                      '%.3f',    1;   ...
            'Accuracy p-value',       acc.accuracy_p,             '%.4g',    1};
    else
        rows = { ...
            'Prediction-outcome r',   acc.prediction_outcome_r,   '%.3f', 1; ...
            'Predicted R^2',          acc.predicted_r2,           '%.3f', 1; ...
            'Out-of-sample R^2',      acc.out_of_sample_r2,       '%.3f', 1; ...
            'RMSE',                   acc.rmse,                   '%.4g', 1; ...
            'MAE',                    acc.mae,                    '%.4g', 1; ...
            'MSE',                    acc.mse,                    '%.4g', 1};
    end
    for i = 1:size(rows, 1)
        v = rows{i, 2};
        if isempty(v) || (isnumeric(v) && isnan(v)), continue; end
        fprintf(['    %-24s ' rows{i, 3} '\n'], rows{i, 1}, rows{i, 4} * v);
    end
    if ~isnan(acc.n)
        fprintf('    %-24s %d\n', 'N', acc.n);
    end
    fprintf('\n');
end


% -------------------------------------------------------------------------
function v = first_num(x)
    if isempty(x), v = NaN; elseif isnumeric(x), v = x(1); else, v = NaN; end
end

function v = fill(v, alt)
    if isnan(v), v = alt; end
end

function out = ternary(tf, a, b)
    if tf, out = a; else, out = b; end
end
