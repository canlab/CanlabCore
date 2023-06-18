function sig = holm_sidak(pVector, alpha)
% performs Holm-Sidak correction for multiple independent comparisons. Sidak correction is a slightly less conservative bonferonni
% Holm's method is a step down approach implemented by the while loop iterations. More commonly used in the context of bonferonni
# rather than Sidak, but equally valid in both cases.
% 
% pVector - array of p-values obtained from multiple independent tests
% alpha - Overall desired false positive rate, e.g. 0.05.
% sig - boolean array of length(pVector) indicating whether a pVector element is significant
%
% Bogdan Petre, 6/19/2023
%
% WARNING: This code is an alpha version, not fully tested for all edge cases

[sortedP, argsort] = sort(pVector, 'descend'); % pValues sorted from largest to smallest, and indices mapping pVector to sortedP

sig = zeros(length(sortedP),1);
while sortedP(end) < 1 - (1-alpha)^(1/length(sortedP))
    sig(argsort(end)) = 1;
    sortedP(end) = [];
    argsort(end) = [];
end

end
