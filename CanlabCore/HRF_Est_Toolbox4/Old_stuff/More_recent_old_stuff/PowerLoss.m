function [PowLoss] = PowerLoss(modres, modfit, moddf, tc, TR, Run, alpha)
% function [PowLoss] = PowerLoss(modres, modfit, moddf, tc, TR, Run, alpha)
%
% Estimates Power-loss due to mis-modeling.
%
% INPUT:
%
% modres - residuals
% modfit - model fit
% moddf  - model degrees of freedom
% tc     - time course
% TR     - time resolution
% Runs   - expermental design
% alpha  - alpha value
%
% OUTPUT:
%
% PowLoss - Estimated power loss
%
%

len = length(tc);               % length of time course
T = round(30./TR);              % length of estimated HRF
tstar = tinv(1-alpha,moddf);    % t-threshold

% Fit FIR model to find 'baseline' power.
[h, fit, e] = Fit_sFIR(tc,TR,Run,T,1);
s = (1/(len-T))*e'*e;
t = 1/sqrt(s*inv(fit'*fit));
basePow = 1- nctcdf(tstar,(len-T),t);

% Compute model power.
sig = (1/moddf)*modres'*modres;
ts = 1/sqrt(sig*inv(modfit'*modfit));
modPow = 1- nctcdf(tstar,moddf,ts);

% Compute 'power loss'
PowLoss = basePow - modPow;
