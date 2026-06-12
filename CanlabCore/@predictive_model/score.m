function v = score(obj, X, Y)
% score  Evaluate obj.scorer on (Y, predict(obj, X)).
%
% :Usage:
% ::
%     v = score(pm, X, Y);
%
% If obj.scorer is empty, picks the default for obj.task via
% cv_scorer.default_for_task (balanced_accuracy for classification,
% r2 for regression).
%
% If the scorer needs continuous scores (e.g. roc_auc), pulls them
% from the second output of predict(obj, X); otherwise just uses
% the discrete predictions.
%
% :Inputs:
%
%   **obj:**
%        a fitted @predictive_model.
%
%   **X:**
%        [n x p] predictor matrix in the model's feature space.
%
%   **Y:**
%        [n x 1] true outcomes to score predictions against.
%
% :Outputs:
%
%   **v:**
%        scalar score from obj.scorer (higher is better unless the scorer's
%        greater_is_better is false, e.g. rmse).
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm');
%     X = dat.dat'; Y = dat.Y;
%     pm = predictive_model('algorithm','svm','task','classification');
%     pm = fit(pm, X, Y);
%     v  = score(pm, X, Y);            % in-sample accuracy (optimistic!)
%
% :See also:
%   fit, predict, crossval, cv_scorer

    if isempty(obj.scorer)
        if isempty(obj.task)
            error('predictive_model:score:NoTask', ...
                'Set obj.task or obj.scorer before calling score().');
        end
        obj.scorer = cv_scorer.default_for_task(obj.task);
    end

    if isobject(obj.scorer) && obj.scorer.needs_continuous
        [yhat, scores] = predict(obj, X);
        v = obj.scorer.score(Y, yhat, scores);
    else
        yhat = predict(obj, X);
        v = obj.scorer.score(Y, yhat);
    end
end
