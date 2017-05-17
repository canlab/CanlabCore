  function [pval]=correst(v1,v2,ntrials,ntrain,ntest,alpha)

% [pval]=correst(v1,v2,ntrials,ntrain,ntest,alpha)
%
% calculates the corrected resampled ttest of Naudeau & Bengio 
% v1 & v2 are vectors of error rates of two algorithms over ntrials trials

  
  quant = tinv(1-alpha/2,ntrials-1);
  
 %[ abs(mean(v1)-mean(v2))  quant*sqrt( (1/ntrials+ntest/ntrain)*(var(v1-v2)) ) ]
  
    if abs(mean(v1)-mean(v2)) > quant*sqrt( (1/ntrials+ntest/ntrain)*(var(v1-v2)) ) 
        pval = 0; 
    else
        pval = 1;
    end