function [yfit, vox_weights, intercept] = cv_mlpcr_bt(xtrain, ytrain, xtest, cv_assignment, varargin)    
    [~,B,~] = mlpcr2(xtrain, ytrain,  varargin{:});
                
    intercept = B(1);
    vox_weights = B(2:end);
    yfit = xtest*vox_weights + intercept;
end