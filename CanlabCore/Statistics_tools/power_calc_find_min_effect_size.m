function [d_max, d_max_2group] = power_calc_find_min_effect_size(n, pthr)
% Find effect size detectable with 80% power for a given sample size (n) for 1-group t-test and 2-group comparisons
% see power_calc.m
% 
% Tor Wager, 2023
%
% % Example: 
% For n = 120, find power to detect effect in a planned test at p < 0.05 two-tailed
% [d_max, d_max_2group] = power_calc_find_min_effect_size(120, 0.025);
%
% % ... and in a standard "guesstimate" of the potential FDR-corrected threshold at q < 0.05, using p < 0.001 one tailed (p < 0.002)
% [d_max, d_max_2group] = power_calc_find_min_effect_size(120, 0.002);

dvals = [.001:.01:2]';
obspow = zeros(size(dvals, 1), 2);

for i = 1:length(dvals)
[~, ~, obspow(i, :)] = power_calc(dvals(i), pthr, n, 'd', false, false, false, 'noverbose');
end
d_max = dvals(find(obspow(:, 1) >= 0.8));
d_max = d_max(1);

d_max_2group = dvals(find(obspow(:, 2) >= 0.8));
d_max_2group = d_max_2group(1);

fprintf(['With n = %3.0f and p < %0.3f, we have 80%% power to detect effects of d = %3.2f or larger in 1 group, and d = %3.2f or larger in a two-group comparison\n'], n, pthr, d_max, d_max_2group);


