function out = permutation_test_site_specificity(RSA, varargin)
% permutation_test_site_specificity  Permutation test for cross-subject bodysite
% specificity that respects the dependency structure of the RSA.
%
% NULL HYPOTHESIS & EXCHANGEABILITY
%   H0: bodysite labels carry no cross-subject representational information.
%   Under H0, the site labels are exchangeable WITHIN each subject (we must keep
%   subject identity fixed, because subjects differ in overall map similarity /
%   SNR). Each permutation independently shuffles the site labels within every
%   subject, then recomputes the SUBJECT-LEVEL specificity statistic
%   (mean over anchor subjects of same-site minus different-site cross-subject
%   similarity, in Fisher-z space for the correlation metric). The observed
%   statistic is compared to this null distribution.
%
%   This is more defensible than ttest2 on pooled pairwise cells (which ignores
%   non-independence) and complements subject_level_crosssubject_effect (which
%   assumes normality of the per-subject contrast).
%
% :Usage:
% ::
%   out = permutation_test_site_specificity(RSA)
%   out = permutation_test_site_specificity(RSA, 'nperm', 5000, 'seed', 1)
%
% :Inputs:
%   **RSA:** struct from build_crosssubject_signature_rsa.
%
% :Optional Inputs:
%   **'nperm':**     number of permutations (default 5000).
%   **'seed':**      RNG seed for reproducibility (default 0). [] = do not set.
%   **'tail':**      'right' (default; same>diff), 'both', 'left'.
%   **'doverbose':** print summary (default true).
%
% :Outputs:
%   **out:** struct with
%      .observed      observed subject-level specificity statistic (z-space)
%      .null          nperm x 1 null statistics
%      .p             permutation p-value (with +1 correction)
%      .z             (observed - mean(null)) / std(null)
%      .nperm, .tail, .metric
%
% :Examples:
% ::
%   out = permutation_test_site_specificity(RSA_dp_corr, 'nperm', 5000);
%   fprintf('specificity p_perm = %.4g\n', out.p);
%
% :See also: subject_level_crosssubject_effect
%
% ..
%    2026 Michael Sun. GPL v3.
% ..

p = inputParser;
p.addParameter('nperm', 5000, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('seed', 0, @(x) isempty(x) || isnumeric(x));
p.addParameter('tail', 'right', @(x) any(strcmpi(x, {'right','both','left'})));
p.addParameter('doverbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
nperm     = round(p.Results.nperm);
seed      = p.Results.seed;
tail      = lower(p.Results.tail);
doverbose = logical(p.Results.doverbose);

if ~isempty(seed), rng(seed); end

nSub   = RSA.nSub;
M      = RSA.M;
iscorr = strcmp(RSA.metric, 'correlation');
if iscorr, Mt = rsa_fisher_z(M); else, Mt = M; end

subj = RSA.subject_idx;
site = RSA.site_idx;

% Observed statistic with true site labels.
observed = local_subject_spec(Mt, subj, site, nSub);

% Permutation: shuffle site labels within each subject.
nullstat = nan(nperm,1);
% Precompute per-subject node indices for fast within-subject shuffling.
subj_nodes = arrayfun(@(s) find(subj == s), 1:nSub, 'UniformOutput', false);

for it = 1:nperm
    site_perm = site;
    for s = 1:nSub
        idx = subj_nodes{s};
        site_perm(idx) = site(idx(randperm(numel(idx))));
    end
    nullstat(it) = local_subject_spec(Mt, subj, site_perm, nSub);
end

% Permutation p-value (+1 correction).
switch tail
    case 'right'
        pval = (1 + sum(nullstat >= observed)) / (nperm + 1);
    case 'left'
        pval = (1 + sum(nullstat <= observed)) / (nperm + 1);
    case 'both'
        pval = (1 + sum(abs(nullstat) >= abs(observed))) / (nperm + 1);
end

out = struct();
out.observed = observed;
out.null     = nullstat;
out.p        = pval;
out.z        = (observed - mean(nullstat, 'omitnan')) / std(nullstat, 'omitnan');
out.nperm    = nperm;
out.tail     = tail;
out.metric   = RSA.metric;

if doverbose
    fprintf('\nPermutation test of cross-subject site specificity (metric=%s)\n', RSA.metric);
    fprintf('  observed same-diff = %.4f | null mean = %.4f (SD %.4f)\n', ...
        observed, mean(nullstat,'omitnan'), std(nullstat,'omitnan'));
    fprintf('  permutation p = %.4g (%s, %d perms), z = %.2f\n', pval, tail, nperm, out.z);
end

end


function stat = local_subject_spec(Mt, subj, site, nSub)
% Subject-level same-site minus different-site, averaged over anchor subjects.
spec = nan(nSub,1);
for a = 1:nSub
    is_a  = subj == a;
    not_a = subj ~= a;
    sites_a = unique(site(is_a));

    same_vals = [];
    diff_vals = [];
    for b = sites_a(:)'
        a_node = is_a & (site == b);
        same_p = not_a & (site == b);
        diff_p = not_a & (site ~= b);

        sv = Mt(a_node, same_p); sv = sv(isfinite(sv));
        dv = Mt(a_node, diff_p); dv = dv(isfinite(dv));
        same_vals = [same_vals; sv(:)]; %#ok<AGROW>
        diff_vals = [diff_vals; dv(:)]; %#ok<AGROW>
    end
    spec(a) = mean(same_vals,'omitnan') - mean(diff_vals,'omitnan');
end
stat = mean(spec, 'omitnan');
end
