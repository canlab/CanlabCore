function res = lasso_selection(lasso_obj, method)
% -------------------------------------------------------------------------
% function res = lasso_selection(lasso_obj, method)
% -------------------------------------------------------------------------
% PURPOSE:
% This function uses Zou and Hastie's unbiased estimate of the 
% degrees of freedom of the lasso to select the amount of regularization
% -------------------------------------------------------------------------
% INPUTS:
% lasso_obj:  a structure as returned by the lasso function
% method: one of the following strings (not case sensitive):
%         "aic": uses Akaike information criterion (AIC)
%         "bic": uses Schwartz's Bayesian information criterion (BIC)
%         "all": returns all selecctions above
% -------------------------------------------------------------------------
% OUTPUTS:
% A copy of the lasso_obj structure with the following fields added 
% according to the choice of method:
% AIC/BIC/gMDL/AICC are structures containing:
%      path:      the path of the criterion statistic
%      minidx:    the index at which the criterion is minimized
%      nbetas:    the normalized coefficients picked by the criterion
%      betas:     the coefficients picked by the criterion
%      intercept: the intercept picked by the criterion
%      lambda:    the intercept picked by the criterion
% -------------------------------------------------------------------------
% Author: Peng Zhao and Guilherme V. Rocha
%         Department of Statistics
%         University of California, Berkeley
%         gvrocha@stat.berkeley.edu, gvrocha@gmail.com
% 2006/09
% -------------------------------------------------------------------------
%
% SEE ALSO lasso

res   = lasso_obj;

y = res.ymean+res.yscale*res.ny;

[junk, residuals] = lasso_predict(res, [], res.normalized_penalty);
tss               = (y-mean(y))'*(y-mean(y));
rss               = sum(residuals.^2);
df                = res.df;
loglikelihood     = 0.5*res.sample_size*log(rss);

if strcmpi(method,'AIC')|strcmpi(method,'all')
  res.AIC.path           = loglikelihood+2*res.dfs;
  [dummy,res.AIC.minidx] = min(res.AIC.path);                       
  res.AIC.beta           = res.beta(res.AIC.minidx,:);                  
  res.AIC.intercept      = res.intercept(res.AIC.minidx,:);             
  res.AIC.nbeta          = res.nbeta(res.AIC.minidx,:);                     
  res.AIC.lambda         = res.lambda(res.AIC.minidx,:);                     
end;

if strcmpi(method,'BIC')|strcmpi(method,'all')
  res.BIC.path           = loglikelihood+log(res.sample_size)*res.dfs;
  [dummy,res.BIC.minidx] = min(res.BIC.path);                          
  res.BIC.beta           = res.beta(res.BIC.minidx,:);                     
  res.BIC.intercept      = res.intercept(res.BIC.minidx,:);               
  res.BIC.nbeta          = res.nbeta(res.BIC.minidx,:);                     
  res.BIC.lambda         = res.lambda(res.BIC.minidx,:);                     
end;

if strcmpi(method,'gMDL')|strcmpi(method,'all')
  warning off;
  F                                  = ((tss-rss)./rss).*((result(iii).sample_size-df)./df);
  warning on;
  res.gMDL.weights        = 0.5*log(F);
  res.gMDL.path           = loglikelihood+res.gMDL.weights.*res.dfs;
  [dummy,res.gMDL.minidx] = min(res.gMDL.path);                          
  res.gMDL.beta           = res.beta(res.gMDL.minidx,:);                     
  res.gMDL.intercept      = res.intercept(res.gMDL.minidx,:);               
  res.gMDL.nbeta          = res.nbeta(res.gMDL.minidx,:);                     
  res.gMDL.lambda         = res.lambda(res.gMDL.minidx,:);                     
end;

if strcmpi(method,'AICC')|strcmpi(method,'all')
  res.AICC.penalty       = 0.5*sample_size*(1+(df/sample_size))./(1-((df+2)/sample_size));
  res.AICC.path           = loglikelihood+res.AICC.penalty;
  [dummy,res.AICC.minidx] = min(res.AICC.path);                          
  res.AICC.beta           = res.beta(res.AICC.minidx,:);                     
  res.AICC.intercept      = res.intercept(res.AICC.minidx,:);               
  res.AICC.nbeta          = res.nbeta(res.AICC.minidx,:);                     
  res.AICC.lambda         = res.lambda(res.AICC.minidx,:);                     
end;
