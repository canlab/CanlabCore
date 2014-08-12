function [h,sig,ci,tval,ser] = t_test2(x,m,alpha,tail)
%TTEST  Hypothesis test: Compares the sample average to a constant.
%   [H,SIG] = TTEST(X,M,ALPHA,TAIL) performs a T-test to determine
%   if a sample from a normal distribution (in X) could have mean M.
%   M = 0, ALPHA = 0.05 and TAIL = 0 by default.
%
%   The Null hypothesis is: "mean is equal to M".
%   For TAIL=0,  alternative: "mean is not M".
%   For TAIL=1,  alternative: "mean is greater than M"
%   For TAIL=-1, alternative: "mean is less than M"
%   TAIL = 0 by default.
%
%   ALPHA is desired significance level. 
%   SIG is the probability of observing the given result by chance
%   given that the null hypothesis is true. Small values of SIG cast
%   doubt on the validity of the null hypothesis.
%   H=0 => "Do not reject null hypothesis at significance level of alpha."
%   H=1 => "Reject null hypothesis at significance level of alpha."

%   References:
%      [1] E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, 1970, page 206. 

%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 2.7 $  $Date: 1998/05/28 20:13:56 $

if nargin < 1, 
    error('Requires at least one input argument.'); 
end

[m1 n1] = size(x);
if (m1 ~= 1 & n1 ~= 1) 
    error('First argument has to be a vector.');
end

if nargin < 2
    m = 0;
end
 
if nargin < 4, 
    tail = 0; 
end 

if nargin < 3, 
    alpha = 0.05; 
end 

if (alpha <= 0 | alpha >= 1)
    fprintf('Warning: significance level must be between 0 and 1\n'); 
    h = NaN;
    sig = NaN;
    ci = [NaN NaN];
    return;
end

samplesize  = length(x);
xmean = mean(x);
ser = std(x) ./ sqrt(samplesize);
tval = (xmean - m) / ser;
sig = tcdf(tval,samplesize - 1);

% the significance just found is for the  tail = -1 test

crit = tinv(1 - alpha,samplesize - 1) .* ser;

if tail == 1
    sig = 1 - sig;
elseif tail == 0
    sig = 2 * min(sig,1 - sig);
    crit = tinv((1 - alpha / 2),samplesize - 1) .* ser;
end

ci = [(xmean - crit) (xmean + crit)];

% if tval > 0, sig = 1-sig;,end

% Determine if the actual significance exceeds the desired significance
h = 0;
if sig <= alpha, 
    h = 1; 
end 

if isnan(sig), 
    h = NaN; 
end
