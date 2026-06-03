function obj = calibrate(obj, X, Y, varargin)
% calibrate  Fit Platt scaling for probabilistic classification outputs.
%
% Platt scaling fits a sigmoid to the model's decision-function scores
% so that the resulting probability is well-calibrated against observed
% class frequencies:
%     P(y = positive | s) = 1 / (1 + exp(-(A*s + B)))
%
% :Usage:
% ::
%     pm = predictive_model('algorithm','svm','task','classification');
%     pm = calibrate(pm, X, Y);                       % cv internally to get unbiased s
%     pm = calibrate(pm, X, Y, 'method', 'isotonic'); % isotonic regression instead
%
%     P = predict_proba(pm, X_new);
%
% Pipeline:
%   1. If obj.fit_type is not already 'crossval', run crossval() so
%      the held-out continuous scores aren't optimistically biased.
%   2. Fit the sigmoid (or isotonic regression) on (cv_scores, Y_binary).
%   3. Store the calibrator in obj.fitted_values.calibrator.
%
% After: predict_proba(pm, X_new) returns calibrated probabilities of
% the positive class.
%
% :Inputs:
%
%   **obj:**
%        a @predictive_model (classification task).
%
%   **X:**
%        [n x p] predictor matrix.
%
%   **Y:**
%        [n x 1] binary outcome.
%
% :Optional Inputs (name/value):
%   'method'   'platt' (default) or 'isotonic'
%
% :Outputs:
%
%   **obj:**
%        the @predictive_model with obj.fitted_values.calibrator populated
%        (Platt A/B parameters or isotonic knots) and obj.do_calibrate set.
%        Use predict_proba(obj, X_new) afterwards.
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm');
%     X = dat.dat'; Y = dat.Y;
%     pm = predictive_model('algorithm','svm','task','classification');
%     pm = calibrate(pm, X, Y, 'method', 'platt');
%     P  = predict_proba(pm, X);       % well-calibrated P(class = +1)
%
% :See also:
%   predict_proba, crossval, predict

    pi = inputParser; pi.KeepUnmatched = true;
    addParameter(pi, 'method', 'platt');
    parse(pi, varargin{:});
    method = lower(pi.Results.method);

    if isempty(obj.task), obj.task = 'classification'; end
    if ~strcmpi(obj.task, 'classification')
        error('predictive_model:calibrate:NotClassification', ...
            'calibrate only valid for classification (got task=''%s'').', obj.task);
    end

    % 1. Get unbiased held-out continuous scores via crossval.
    if ~strcmp(obj.fit_type, 'crossval') || isempty(obj.fitted_values) ...
            || ~isfield(obj.fitted_values, 'scores') ...
            || isempty(obj.fitted_values.scores)
        obj = crossval(obj, X, Y);
    end

    scores  = obj.fitted_values.scores(:);
    Y_use   = obj.Y(:);
    classes = unique(Y_use);
    if numel(classes) ~= 2
        error('predictive_model:calibrate:Binary', ...
            'calibrate currently supports binary classification only.');
    end
    pos = max(classes);
    y_bin = double(Y_use == pos);

    valid = ~isnan(scores) & ~isnan(y_bin);
    s = scores(valid); yb = y_bin(valid);

    switch method
        case 'platt'
            % Logistic regression: y_bin ~ sigmoid(A*s + B).
            % glmfit returns [intercept; slope] for the canonical link.
            b = glmfit(s, yb, 'binomial', 'link', 'logit');
            B = b(1);
            A = b(2);
            calibrator = struct('method', 'platt', 'A', A, 'B', B, ...
                                'positive_class', pos);

        case 'isotonic'
            % Non-parametric monotone fit via pool-adjacent-violators.
            [s_sorted, ord] = sort(s);
            y_sorted = yb(ord);
            iso = predictive_model.pava(y_sorted);
            calibrator = struct('method', 'isotonic', ...
                                's_knots', s_sorted, ...
                                'p_knots', iso, ...
                                'positive_class', pos);
        otherwise
            error('predictive_model:calibrate:UnknownMethod', ...
                'Unknown method ''%s''. Options: ''platt'', ''isotonic''.', method);
    end

    obj.fitted_values.calibrator = calibrator;
    obj.do_calibrate = true;

    if strcmp(method, 'platt')
        obj.history{end+1, 1} = sprintf( ...
            'calibrate (Platt): A=%.3f, B=%.3f', calibrator.A, calibrator.B);
    else
        obj.history{end+1, 1} = sprintf( ...
            'calibrate (isotonic): %d monotone knots', numel(calibrator.p_knots));
    end
end
