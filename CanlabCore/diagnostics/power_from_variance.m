function [power, t, d] = power_from_variance(con, N, sig2b, sig2wi, pthresh)
% Power and effect size measures, given contrast, N, and variance component
% estimates
%
% :Inputs:
%
%   **con:**
%        contrast/effect magnitude estimate; "mean difference"
%
%   **N:**
%        sample size
%
%   **sig2b:**
%        between-subjects variance estimate
%
%   **sig2wi:**
%        within-subjects variance estimate
%           *note* this is not the "raw" within-subjects variance; it is
%           the contribution to the group (2nd-level) variance, which is
%           sig2within / number of images within-subjects
%
%   **pthresh:**
%        alpha (Type I error) rate; p-value threshold for power
%             calculation
%
% con, sig2b, and sig2wi can all be vectors, so you can run this function
% voxel-wise for a whole map at once
%
% :Outputs:
%
%   **power:**
%        Power from 0 to 1
%
%   **t:**
%        effect size : expected t-value
%
%   **d:**
%        effect size : Cohen's d
%
% See canlab_effect_size_map for a whole-brain, image-based power mapping
% function, and canlab_power_allocation_curves for normative power curves.
%
% Note: power is computed with the classical "shifted central t"
% approximation, 1 - tcdf(u - t, N - 1). canlab_effect_size_map uses the
% noncentral t distribution by default, which is exact under the model.
%
% ..
%    Tor Wager, August 2009
%    2026-09: Critical t value now uses N - 1 degrees of freedom (was N),
%    matching the one-sample t-test df used for the power calculation.
% ..

% t-value threshold for significance at alpha level pthresh (two-tailed),
% for a one-sample t-test with N - 1 degrees of freedom
u_unc = tinv(1 - pthresh./2, N - 1);

% expected t-value based on N, con, variance components
t = ( abs(con) .* sqrt(N) ) ./ (sig2b + sig2wi) .^ .5;

d = t ./ sqrt(N);

% power at this alpha level
power = 1 - tcdf(u_unc - t, N - 1);


end
