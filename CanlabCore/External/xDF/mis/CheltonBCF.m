function [V,Z,P,BCF] = CheltonBCF(Y,T,CF,verbose)
% Estimating and correcting the degrees of freedom based on approximations
% presented in Barteltt's 1946 and Quenouille's 1947.
%
% Pyper calls this Chelton's method, but I can't find the paper for the
% life of me and the equation is identical to the Bartlett's. Also,
% Bartlett's paper is much older than Chelton's.
%
%%%INPUTS
%   Y : Time series.
%   T : #data points, just for sanity check
%   CF: Curbing Factor
%
%%%OUTPUTS
%   V : Variance
%   Z : Z-score (via Fisher's transformation)
%   P : p-value (two-sided on a-level: 5%)
%   BCF : Correction Factor for effective degrees of freedom
%
%%%REFERENCE
%   Bartlett, M. S. (1946). On the Theoretical Specification and Sampling 
%   Properties of Autocorrelated Time-Series. Supplement to the Journal of 
%   the Royal Statistical Society, 8(1), 27. http://doi.org/10.2307/2983611
%   
%   Quenouille, M. H. (1947). Notes on the Calculation of Autocorrelations 
%   of Linear Autoregressive Schemes. Biometrika, 34(3/4), 365. 
%   http://doi.org/10.2307/2332450
%
%   Bayley, G. V., & Hammersley, J. M. (1946). The ?Effective? Number of 
%   Independent Observations in an Autocorrelated Time Series. Supplement 
%   to the Journal of the Royal Statistical Society, 8(2), 184. 
%   http://doi.org/10.2307/2983560
%
%   Pyper, B. J., & Peterman, R. M. (1998). Comparison of methods to 
%   account for autocorrelation in correlation analyses of fish data, 
%   2140, 2127?2140.
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2018
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________
if ~exist('verbose','var'); verbose = 1;  end

if ~exist('CF','var') || isempty(CF)
    if verbose; disp('CF hasnt been assigned, so set to 1.'); end;
    CF = 1; 
end

if size(Y,1)~=T; Y=Y'; end; %IxT


ac        = AC_fft(Y,T);
ac(:,[1 end]) = []; %remove the 0lag, we'll take care of it later.

ac = curbtaperme(ac,T-2,round((T-2)./CF)); %curb according to Pyper and Peterman

BCF       = 1+2.*(ac*ac');

if any(any(BCF>T)); BCF (BCF>T) = T; end
if any(any(BCF<1)) && verbose; disp('BCF was below 1, set back to 1.'); BCF(BCF<1) = 1; end;

V         = BCF./T; %is it?
Z         = atanh(corr(Y)).*sqrt(T./BCF-3);
P         = 2 * normcdf(-abs(Z));  

end

function ct_ts=curbtaperme(ts,T,M)
% Curb the autocorrelations, according to Anderson 1984
% multi-dimensional, and therefore is fine!
%SA, Ox, 2018
    if ~sum(ismember(size(ts),T)); error('There is something wrong, mate!'); end
    if size(ts,2) ~= T; ts = ts'; end
    
    M          = round(M);
    msk        = zeros(size(ts));
    msk(:,1:M) = 1;
    ct_ts      = msk.*ts;
end