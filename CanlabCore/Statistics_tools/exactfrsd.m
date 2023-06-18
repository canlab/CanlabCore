% Nonparametric test for exact likelihoods of pairwise contrasts in a 
% repeated measures design. Is to the Friedman test as a Tukey test is 
% to an rm-ANOVA.
% function p = exactfrsd(d,k,n)
%
% Output ::
%
%   p 	- Likelihood under the null (both groups equal). aka pValue
%
% Input ::
%
%   d 	- sum of Friedman rank difference. These can be obtained with 
%       	a = multcompare(friedman(data)); ranksum = n*a(:,4)
%   k	- number of treatments (e.g. classifiers) tested
%   n 	- number of groups (e.g. datasets) on which repeated measures were 
%           evaluated
%
% adapted from R code. Cite 
% Eisinga R, Heskes T, Pelzer B, Grontenhuis MT (2017) "Exact p-values for 
%   pairwise comparison of Friedman rank sums, with applications to 
%   comparing classifiers". BMC Bioinformatics 18:68.
%
% shamelessly copied from the author's published R script by Bogdan Petre
% April, 2020
%
% Warning: These computations can involve large numbers. The R script (and
% the original paper) solves numerical overflow problems using an R package
% for arbitrary precision arithmetic. We use matlab's vpa() function but
% this is not tested carefully. If in doubt I would defer to the R code.
function p = exactfrsd(d, k, n)
    if n < 1, error('n >= 1 please'); end
    if k < 2, error('n >= 2 please'); end
    if d > n*(k-1), error('|d| <= n*(k-1) please'); end
    
    assert(length(d) == 1, 'Please run on scalar value d');
    
    if d < 0, d=-1*d; end
    
    % rscript expresses variable precision arithmetic bounds in bits, but
    % matlab's vpa() function expects digits. Lets do the conversion for
    % equivalent performance
    bits2048 = ceil(2048*log(2)/log(10));
    
    result = 0;
    for h=0:n
        sum1 = nchoosek(vpa(n, bits2048),h) / ( vpa(k, bits2048)^h * vpa((1-k), bits2048)^n );
        sum2 = 0;
        for s=0:h
            if any(k*s-d-h >= 0)
                sum2 = sum2 + (-1)^s * nchoosek(vpa(2*h, bits2048),h+s) * nchoosek(vpa(k*s-d+h, bits2048), k*s-d-h);
            end
        end
        result = result + sum1*sum2;
    end
        
    if d == 0
        p = 1;
    else
        % vpa() no longer needed, so convert back to double
        p = double(2*result);
    end
end
