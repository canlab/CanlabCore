function [AICvalue] = AIC_all(modRes, num_vox, num_time, k)
%function [PowLoss] = PowerLoss_all(modRes, modFit, moddf, TC, TR, Runc, alpha, Z)


% function [PowLoss] = PowerLoss(modres, modfit, moddf, tc, TR, Run, alpha)
%
% Estimates Power-loss due to mis-modeling.
%
% INPUT:
%
% modRes - residuals (voxels x time course)
% modFit - model fit (voxels x time course)
% moddf  - model degrees of freedom
% fullRes - residuals for full model (voxels x time course)
% fullFit - model fit for full model (voxels x time course)
% fulldf  - model degrees of freedom for full model 
% TR      - time resolution
% num_vox - number of voxels
% num_time - number of volumes
% alpha  - alpha value
%
% OUTPUT:
%
% PowLoss - Estimated power loss
%
%


AICvalue = zeros(num_vox,1);



% Fit FIR model to find 'baseline' power.
%[h, fit, e] = Fit_sFIR(tc,TR,Run,T,1);
%[~, Fit, Resid, ~, ~] = Fit_sFIR_all(TC, TR, Runc, 30, 0, Z);


for v=1:num_vox

    mres = modRes(v,:)';
    s = (1/(num_time - k))*(mres')*mres;
    AICvalue(v) = 2*k + num_time*log(2*pi) + num_time*log(s^2) + (1/s^2)*(mres')*mres;

end
