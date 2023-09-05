function [p sres sres_ns] = ResidScan(res, FWHM)
% function [p sres sres_ns] = ResidScan(res, FWHM)
%
% Calculates P(M>=t) where M is the max value of the smoothed residuals.
% In this implementation the residuals are smoothed using a Gaussian
% kernel.
% 
% INPUT:
%
% res - residual time course
% FWHM - Full Width Half Maximum (in time units)
%
% OUTPUT:
%
% p - pvalues
% sres - smoothed residuals
% sres_ns - smoothed residuals (non standardized) 
%
% By Martin Lindquist & Ji-Meng Loh, July 2007
%
% Edited by ML on 10/02/09

res_ns = res;
res = res./std(res);
len = length(res);

% Create Gaussian Kernel
sig = ceil(FWHM/(2*sqrt(2*log(2))));    
klen = 3*sig;
kern = normpdf((-klen:klen),0,sig); 
kern = kern./sqrt(sum(kern.^2));

% Convolve
x = conv(res,kern);
sres = x((klen + 1):(end-klen));

x = conv(res_ns,kern/sum(kern));
sres_ns = x((klen + 1):(end-klen));


% Find Max value
[a,location] = max(abs(sres));

% Find p-values using Gaussian Random Field theory
z = Euler_p(1, a, len, FWHM);
z = 2*z;        %Two-sided test
p = min(1, z);

end

% END MAIN FUNCTION


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subfunctions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function pval = Euler_p(myDim, value, N, fwhm)
% function z = Euler_p(myDim, value, N, fwhm)
%
% Finds the p value using the expected Euler characteristic. 
% 
% This function returns P(M \ge value) using the approximation 
% \sum_{d=0}^D R_d(V) \rho_d(value) following Worsley et al's "A Unified
% Statistical Approach for Determining Significant Signals in Images of
% Cerebral Activation".
%
% INPUTS:
%
% myDim - the number of dimensions in the data
% value - the value of the maximum. 
% N     - the number of (time) points in that 1 dimension 
% fwhm  - the full width half maximum
%
% OUTPUTS:
%
% pval - the p-value

% NOTE: CURRENTLY THIS FUNCTION IS ONLY IMPLEMENTED FOR THE 1D CASE

  % Constants 
  myfactor = 4*log(2);
  pi2 = 2*pi;
  exptsq = exp(-(value^2)/2);

  % Euler Characteristc Densties
  rho = zeros(5,1);
  rho(1) = 1-normcdf(value);
  rho(2) = myfactor^(0.5)*exptsq/pi2;
  rho(3) = myfactor * exptsq * value / (pi2 ^ (1.5));
  rho(4) = myfactor ^ (1.5) * exptsq * (value^2-1) / (pi2 ^2);
  rho(5) = myfactor ^2 * exptsq * (value^3-3*value) / (pi2 ^ (5/2));
     
  % Resel Count
  R0 = 1;
  R1 = N/fwhm;

  % P-value
  pval = R0 * rho(1) + R1 * rho(2);

end


