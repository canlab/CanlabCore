function [yfit, vox_weights, intercept] = cv_mlpcr_wi(xtrain, ytrain, xtest, cv_assignment, varargin)    
    [B0,~,B] = mlpcr2(xtrain, ytrain,  varargin{:});
                
    intercept = B0(1);
    vox_weights = B(2:end);
    yfit = xtest*vox_weights + intercept;
end
