function [yhat, score_out] = predict(obj, X)
% predict  Apply the fitted model to new data.
%
% :Usage:
% ::
%     [yhat, scores] = predict(pm, X_new);
%
% Pipeline:
%   1. Apply the omitted_features mask if any were dropped at fit time.
%   2. Apply the standardize_mu / standardize_sd transform (if fit()
%      standardized X) so the new data lives in the same feature space.
%   3. Delegate to predict(obj.ml_model, X). Second output (continuous
%      scores) is the MATLAB model's native score: signed distance for
%      SVM, posterior probability for fitPosterior-wrapped SVM,
%      predicted value for regressors.

    if isempty(obj.ml_model)
        error('predictive_model:predict:NotFitted', ...
            'Model has not been fitted. Call fit() first.');
    end

    % Apply omitted_features mask if it was set at fit time.
    if islogical(obj.omitted_features) && any(obj.omitted_features)
        X = X(:, ~obj.omitted_features);
    end

    % Apply standardization if fit() recorded one.
    if isstruct(obj.inputParameters) ...
            && isfield(obj.inputParameters, 'standardize_mu')
        X = (X - obj.inputParameters.standardize_mu) ./ ...
             obj.inputParameters.standardize_sd;
    end

    % Delegate to the trained MATLAB model.
    if nargout > 1
        [yhat, score_out] = predict(obj.ml_model, X); %#ok<MINV> external method
    else
        yhat = predict(obj.ml_model, X);
    end
end
