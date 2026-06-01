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
