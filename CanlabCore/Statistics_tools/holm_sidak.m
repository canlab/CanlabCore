function [sig, pthr] = holm_sidak(pVector, alpha)
% performs Holm-Sidak correction for multiple independent comparisons. Sidak correction is a slightly less conservative bonferonni
% Holm's method is a step down approach implemented by the while loop iterations. More commonly used in the context of bonferonni
% rather than Sidak, but equally valid in both cases.
% 
% pVector - array of p-values obtained from multiple independent tests
% alpha - Overall desired false positive rate, e.g. 0.05.
% sig - boolean array of length(pVector) indicating whether a pVector element is significant
%
% Bogdan Petre, 6/19/2023
%
% WARNING: This code is an alpha version, not fully tested for all edge cases
%
% 1.The P values are ranked from smallest to largest.
% 2.Set a value for the significance level, alpha. This is often set to 5%.
% 3.Define k equal to the number of comparisons (length(sortedP) in the code)
% 4.Start with the smallest P value and set i=k. Ask: Is the smallest P value less than the Šídák corrected value 1-(1-alpha)(1/i)?
% 
%   If No: Conclude that none of the comparisons are statistically significant, and you are done.
%   If Yes: Conclude that this comparison is statistically significant, and continue.
% 
% 5.The second to smallest P value is compared next. Set i=K-1. Is the P value less than 1-(1-alpha)(1/i)?
% 
% If No: Conclude that this  comparison (and all with larger P values) is not statistically significant. Exit.
% 
% If Yes: Conclude that this comparison is statistically significant, and continue step 5 with i = K-2, etc.


[sortedP, argsort] = sort(pVector, 'descend'); % pValues sorted from largest to smallest, and indices mapping pVector to sortedP

sig = zeros(length(sortedP),1);

pthr = 1 - (1-alpha)^(1/length(sortedP)); % max P-value for sig results. Anything below this is signficant.

while sortedP(end) < 1 - (1-alpha)^(1/length(sortedP))

    sig(argsort(end)) = 1;
    sortedP(end) = [];
    argsort(end) = [];

    pthr = 1 - (1-alpha)^(1/length(sortedP)); % step-down adjustment of P-value.

end

end
