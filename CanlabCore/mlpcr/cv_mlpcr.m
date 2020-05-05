function [yfit, vox_weights, intercept] = cv_mlpcr(xtrain, ytrain, xtest, cv_assignment, varargin)
    B = mlpcr2open (xtrain, ytrain,  varargin{:});
                
    intercept = B(1);
    vox_weights = B(2:end);
    yfit = xtest*vox_weights + intercept;
end