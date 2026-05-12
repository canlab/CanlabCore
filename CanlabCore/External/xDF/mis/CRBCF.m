function [V,Z,P,BCF] = CRBCF(Y,T,CF,verbose)
% The way that Pyper says is right!
% it is fast function i.e. more than 2 time series are accptable for input 
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
%   Clifford, P., Richardson, S., & Hemon, D. (2016). 
%   Assessing the Significance of the Correlation between Two Spatial 
%   Processes Author ( s ): Peter Clifford , Sylvia Richardson and Denis 
%   Hemon Published by : International Biometric Society Stable 
%   URL : http://www.jstor.org/stable/2532039 
%   REFERENCES Linked , 45(1), 123?134.
%
%   Bayley, G. V., & Hammersley, J. M. (1946). The ?Effective? Number of 
%   Independent Observations in an Autocorrelated Time Series. Supplement 
%   to the Journal of the Royal Statistical Society, 8(2), 184. 
%   http://doi.org/10.2307/2983560
%
%   Pyper, B. J., & Peterman, R. M. (1998). Comparison of methods to 
%   account for autocorrelation in correlation analyses of fish data, 
%   2140, 2127?2140.
%
%  Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols. 
%  "Effective Degrees of Freedom of the Pearson's Correlation Coefficient 
%  under Serial Correlation." bioRxiv (2018): 453795.
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2018
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________
if ~exist('verbose','var')
    verbose = 1; 
end

if ~exist('CF','var') || isempty(CF)
    if verbose; disp('CF hasnt been assigned, so set to 1.'); end;
    CF = 1; 
end

if size(Y,2)~=T
    Y=Y'; %TxI
end

N = size(Y,1);

ac        = AC_fft(Y,T);
ac(:,[1 end]) = []; %remove the 0lag, we'll take care of it later.

ac = curbtaperme(ac,T-2,round((T-2)./CF)); %curb according to Pyper and Peterman

nLg    = T-1;
wgt    = (nLg:-1:2)./(T);

BCF    = 1+2.*(wgt.*ac*ac');

if any(any(BCF>T)); BCF (BCF>T) = T; end
if any(any(BCF<1)) && verbose; disp([num2str(numel(find(BCF<1))) ' BCF was below 1, set back to 1.']); BCF(BCF<1) = 1; end

%atanh(corr(Y)).*sqrt(T./BCF-3)

V            = BCF./T; %is it?
Z            = atanh(corr(Y')).*sqrt(T./BCF-3);
Z(1:N+1:end) = 0;

P            = 2 * normcdf(-abs(Z));
P(1:N+1:end) = 1; 
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