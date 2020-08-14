% computes bootstrap bias corrected accelerated bootstrap confidence
% interval and returns if it includes 0 as a logical vector sig
% mu_est - estimated parameter from complete sample
% bs_est - bootstrap estimates
% jk_est - jacknife estimates
% alpha - desired FPR
%
% Written by Bogdan Petre, 2020
function [sig, CI] = BCa(mu_est, bs_est, jk_est, alpha)
    fprintf('Estimating BCa CIs...\n')
    bias = sum(bs_est < mu_est)./length(bs_est);
    z0 = icdf('norm', bias, 0, 1);

    jk_mean = mean(jk_est);
    num = sum((jk_mean - jk_est).^3);
    den = sum((jk_mean - jk_est).^2);
    a = num ./ (6*den.^(3/2));

    zL = z0 + icdf('norm',alpha/2,0,1);
    alpha1 = normcdf(z0 + zL./(1-a.*zL));
    zU = z0 + icdf('norm',1-alpha/2,0,1);
    alpha2 = normcdf(z0 + zU./(1-a.*zL));
    
    CI = quantile(bs_est, [alpha1, alpha2]);

    sig = prod(sign(CI),2) > 0;
end
