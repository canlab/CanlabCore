function [Z,BCF,P,BCF_tmp] = FoxBCF(Y,T)
% The way that Fox and VanDjk estimate the BCF; i.e. Integral [sci] of ACFs 
% Fox et al 2005 and VanDijek 200X
%
% NB! This eats all the time series, gives you back one universal number!
%
%%%INPUTS
%   Y : Time series.
%   T : #data points, just for sanity check
%
%%%OUTPUTS
%   Z : Z-score (via Fisher's transformation)
%   BCF : Correction Factor for effective degrees of freedom
%
%%%REFERENCES
%   Bartlett, M. S. (1946). On the Theoretical Specification and Sampling 
%   Properties of Autocorrelated Time-Series. Supplement to the Journal of 
%   the Royal Statistical Society, 8(1), 27. http://doi.org/10.2307/2983611
%
%   Dijk, K. R. A. Van, Hedden, T., Venkataraman, A., Evans, K. C., Lazar, 
%   S. W., & Buckner, R. L. (2010). Intrinsic Functional Connectivity As a 
%   Tool For Human Connectomics : Theory , Properties , and Optimization, 
%   2138, 297?321. http://doi.org/10.1152/jn.00783.2009.
%
%   Fox, M. D., Snyder, A. Z., Vincent, J. L., Corbetta, M., Van Essen, D. 
%   C., & Raichle, M. E. (2005). The human brain is intrinsically organized 
%   into dynamic, anticorrelated functional networks. Proceedings of the 
%   National Academy of Sciences of the United States of America, 102(27), 
%   9673?8. http://doi.org/10.1073/pnas.0504136102
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2018
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

if size(Y,1)~=T
    Y=Y'; %IxT
end

ac              = AC_fft(Y,T);
ac(:,[1 end])   = []; %the last one is shit, the first one is 1!!
BCF_tmp         = 1+2.*sum(ac.^2,2);
BCF             = mean(BCF_tmp);
Z               = atanh(corr(Y)).*sqrt(T./BCF-3);
P               = 2 * normcdf(-abs(Z)); 