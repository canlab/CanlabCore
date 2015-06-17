function X = orthogonalize(mX,X,varargin)
% function X = orthogonalize(mX,X,[scale])
% orthogonalizes X with respect to mX, optionally scaling predictors of X% For each nuisance covariate (column of X)
% Regresses out model fits and saves residuals in X


px = mX * pinv(mX);

for i = 1:size(X,2)
    X(:,i) = X(:,i) - px * X(:,i);                      % residuals
    
    if length(varargin) > 0 & std(X(:,i) ~= 0)
        X(:,i) = (X(:,i) - mean(X(:,i))) ./ std(X(:,i));    % re-scale
    end
    
end

return