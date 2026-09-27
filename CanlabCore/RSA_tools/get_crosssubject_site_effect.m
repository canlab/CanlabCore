function [same_z, diff_z, same_r, diff_r] = get_crosssubject_site_effect(RSA)
% get_crosssubject_site_effect  Cross-subject same-site vs different-site cells.
%
% Collects, from the full signature RSA, every CROSS-SUBJECT pair (subject i ~=
% subject j) and splits it into same-bodysite and different-bodysite pools.
% Same-subject and self (diagonal) comparisons are excluded.
%
% IMPORTANT — THIS IS DESCRIPTIVE ONLY.
%   The returned pools are pairwise cells that are NOT independent: each
%   subject's map participates in many pairs, so a two-sample test on these
%   pools (e.g. ttest2(same_z, diff_z)) has a massively inflated effective N and
%   anticonservative p-values. For inference use subject_level_crosssubject_effect
%   (subject-level aggregation) or permutation_test_site_specificity. Use these
%   pools only for descriptive means / bar plots.
%
% Fisher-z (same_z, diff_z) is meaningful ONLY when RSA.metric == 'correlation'.
% For 'cosine'/'dotproduct' the z outputs are returned empty (with a warning),
% since atanh is not defined for unbounded similarities.
%
% :Usage:
% ::
%   [same_z, diff_z, same_r, diff_r] = get_crosssubject_site_effect(RSA)
%
% :Inputs:
%   **RSA:** struct from build_crosssubject_signature_rsa.
%
% :Outputs:
%   **same_r / diff_r:** column vectors of cross-subject same-/different-site
%                        similarities (r-space, for descriptives/plots).
%   **same_z / diff_z:** Fisher-z of the above ('correlation' metric only;
%                        otherwise []).
%
% :See also: subject_level_crosssubject_effect, permutation_test_site_specificity
%
% ..
%    2026 Michael Sun. GPL v3.
% ..

M  = RSA.M;
% Upper triangle, cross-subject only.
ut = triu(true(size(M)), 1);
xs = RSA.cross_subject_mask & ut;

same_pairs = xs &  RSA.same_site_mask;
diff_pairs = xs & ~RSA.same_site_mask;

same_r = M(same_pairs);
diff_r = M(diff_pairs);
same_r = same_r(isfinite(same_r));
diff_r = diff_r(isfinite(diff_r));

if strcmp(RSA.metric, 'correlation')
    same_z = rsa_fisher_z(same_r);
    diff_z = rsa_fisher_z(diff_r);
else
    warning('get_crosssubject_site_effect:metric', ...
        ['Fisher-z is only valid for correlation metric (RSA.metric = ''%s''). ', ...
         'Returning empty z outputs; use r-space descriptives.'], RSA.metric);
    same_z = [];
    diff_z = [];
end

end
