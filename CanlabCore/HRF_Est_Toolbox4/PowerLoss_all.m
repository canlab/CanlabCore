function [PowLoss] = PowerLoss_all(modRes, modFit, moddf, fullRes, fullFit, fulldf, TR, num_vox, num_time, alpha)
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

T = round(30./TR);              % length of estimated HRF
%T = 30;
%num_vox = size(TC,1);
%num_time = size(TC,2);
PowLoss = zeros(num_vox,1);

tstar = tinv(1-alpha,moddf);    % t-threshold

% Fit FIR model to find 'baseline' power.
%[h, fit, e] = Fit_sFIR(tc,TR,Run,T,1);
%[~, Fit, Resid, ~, ~] = Fit_sFIR_all(TC, TR, Runc, 30, 0, Z);


for v=1:num_vox

    fres = fullRes(v,:)';
    ffit = fullFit(v,:)';

    s = (1/(num_time - T))*(fres')*fres;
    t = 1/sqrt(s*inv(ffit'*ffit));
    basePow = 1- nctcdf(tstar,fulldf,t);
    
    % Compute model power.
  
    mres = modRes(v,:)';
    mfit = modFit(v,:)';

    sig = (1/moddf)*(mres')*mres;
    ts = 1/sqrt(sig*inv(mfit'*mfit));
    modPow = 1- nctcdf(tstar,moddf,ts);
    
    % Compute 'power loss'
    PowLoss(v) = basePow - modPow;

end
