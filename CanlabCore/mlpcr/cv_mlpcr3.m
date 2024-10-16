function [yfit, vox_weights, vox_weights_bt, vox_weights_wi, intercept, pc_b, sc_b, pc_w, sc_w, b] = ...
cv_mlpcr3(xtrain, ytrain, xtest, cv_assignment, varargin)

    [B, Bb, Bw, pc_b, sc_b, pc_w, sc_w, b] = mlpcr3(xtrain, ytrain,  varargin{:});
                
    intercept = B(1);
    vox_weights = B(2:end);
    vox_weights_bt = Bb(2:end);
    vox_weights_wi = Bw(2:end);
    try
        yfit = xtest*vox_weights + intercept;
    catch
        keyboard
    end
end