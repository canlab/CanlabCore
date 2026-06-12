function ROC = rocplot(obj, varargin)
% rocplot  ROC curve from a classifier's cross-validated scores.
%
% Thin wrapper around CANlab's roc_plot: feeds the continuous
% cross-validated model scores (obj.fitted_values.dist_from_hyperplane_xval,
% or .scores) and the binarized outcome (obj.Y) into roc_plot.
%
% :Usage:
% ::
%     ROC = rocplot(pm);
%     ROC = rocplot(pm, 'threshold', 0);
%     ROC = rocplot(pm, 'twochoice');      % forced-choice / paired AUC
%
% :Inputs:
%
%   **obj:**
%        a cross-validated @predictive_model (binary classification).
%
% :Optional Inputs (forwarded to roc_plot):
%
%   Any roc_plot option, e.g. 'threshold', 'color', 'twochoice',
%   'include', 'nofig'. See `help roc_plot`. By default a single-interval
%   ROC at the model's natural decision boundary is drawn.
%
% :Outputs:
%
%   **ROC:**
%        the struct returned by roc_plot (AUC, accuracy, sensitivity,
%        specificity, threshold, ...).
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm');
%     pm  = predictive_model('algorithm','svm','task','classification');
%     pm  = crossval(pm, dat.dat', dat.Y);
%     ROC = rocplot(pm);
%     fprintf('AUC = %.3f\n', ROC.AUC);
%
% :See also:
%   roc_plot, plot, confusionchart, report_accuracy, summary, predictive_model

    % Continuous cross-validated scores.
    scores = local_fitted(obj, 'dist_from_hyperplane_xval');
    if isempty(scores), scores = local_fitted(obj, 'scores'); end
    if isempty(scores)
        error('predictive_model:rocplot:NoScores', ...
            ['Need continuous model scores (fitted_values.dist_from_hyperplane_xval ' ...
             'or .scores). Run crossval() on a classification model first.']);
    end
    if isempty(obj.Y)
        error('predictive_model:rocplot:NoY', 'obj.Y is empty.');
    end

    Y = obj.Y(:);
    classes = unique(Y(~isnan(Y)));
    if numel(classes) ~= 2
        error('predictive_model:rocplot:NotBinary', ...
            'rocplot requires a binary outcome (got %d classes).', numel(classes));
    end

    % Binarize: positive = the larger class label (e.g. +1).
    binary_outcome = Y == max(classes);

    valid = ~isnan(scores(:)) & ~isnan(Y);
    ROC = roc_plot(scores(valid), binary_outcome(valid), varargin{:});
end


% -------------------------------------------------------------------------
function v = local_fitted(obj, name)
    if isstruct(obj.fitted_values) && isfield(obj.fitted_values, name)
        v = obj.fitted_values.(name);
    else
        v = [];
    end
end
