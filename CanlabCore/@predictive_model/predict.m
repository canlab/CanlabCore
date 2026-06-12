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
%
% :Inputs:
%
%   **obj:**
%        a fitted @predictive_model (obj.ml_model populated by fit()).
%
%   **X:**
%        [n x p] predictor matrix in the ORIGINAL (pre-omitted-features)
%        feature space; predict() applies the stored omitted-features mask
%        and standardization itself.
%
% :Outputs:
%
%   **yhat:**
%        [n x 1] predicted class labels (classification) or values
%        (regression).
%
%   **score_out:**
%        [n x k] continuous scores: signed distance from the hyperplane
%        (SVM), posterior probability (calibrated SVM), or predicted value
%        (regression). The last column is the positive-class / regression
%        score consumed by crossval and the scorers.
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm');
%     X = dat.dat'; Y = dat.Y;
%     pm = predictive_model('algorithm','svm','task','classification');
%     pm = fit(pm, X, Y);
%     [yhat, scores] = predict(pm, X);
%
% :See also:
%   fit, score, crossval, predict_proba

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

    % PCR (and other struct-backed linear models) predict directly from
    % the stored voxel-space weights: yhat = intercept + X * w.
    if isstruct(obj.ml_model) && isfield(obj.ml_model, 'type') ...
            && strcmp(obj.ml_model.type, 'pcr')
        yhat = obj.ml_model.intercept + X * obj.ml_model.w;
        if nargout > 1, score_out = yhat; end
        return
    end

    % Delegate to the trained MATLAB model.
    % Classification models return [labels, scores]; regression models
    % return just yhat. For regression we treat yhat itself as the
    % "continuous score".
    if nargout > 1
        try
            [yhat, score_out] = predict(obj.ml_model, X); %#ok<MINV> external method
        catch
            yhat = predict(obj.ml_model, X);
            score_out = yhat;
        end
    else
        yhat = predict(obj.ml_model, X);
    end
end
