function out = subject_level_crosssubject_effect(RSA, varargin)
% subject_level_crosssubject_effect  Defensible subject-level test of cross-subject
% same-site > different-site signature similarity.
%
% THE NON-INDEPENDENCE PROBLEM THIS SOLVES
%   The full RSA has nSub*nSites nodes; the cross-subject same-/different-site
%   pools contain hundreds of pairwise cells, but each subject's maps appear in
%   many cells. Treating those cells as independent (ttest2 on the pooled cells)
%   inflates the effective N and the t-statistic dramatically. Here we collapse
%   to ONE observation per subject before testing, so the test has df = nSub-1.
%
% METHOD (leave-one-anchor aggregation)
%   For each anchor subject a:
%     same(a) = mean over all (site b, other-subject) cross-subject pairs of the
%               SAME-site similarity, anchored on a;
%     diff(a) = mean over all cross-subject DIFFERENT-site pairs anchored on a.
%   For the 'correlation' metric these means are taken in Fisher-z space.
%   The contrast same(a) - diff(a) is a single number per subject. A paired
%   (one-sample on the difference) t-test across the nSub anchors tests
%   same > different. Also returns per-anchor, per-site same-site means so two
%   ROIs (e.g. dpIns vs S1) can be compared with a paired t-test across subjects.
%
% :Usage:
% ::
%   out = subject_level_crosssubject_effect(RSA)
%   out = subject_level_crosssubject_effect(RSA, 'tail', 'right')
%
% :Inputs:
%   **RSA:** struct from build_crosssubject_signature_rsa.
%
% :Optional Inputs:
%   **'tail':**      't-test tail for same>diff: 'right' (default), 'both', 'left'.
%   **'doverbose':** print a summary table (default true).
%
% :Outputs:
%   **out:** struct with fields
%      .same_z, .diff_z      nSub x 1 per-anchor means (z-space; raw for non-corr)
%      .spec_z               nSub x 1 same - diff per anchor (the test variable)
%      .same_r, .diff_r      per-anchor means mapped back to r-space
%      .sitewise_same_z      nSub x nSites per-anchor, per-site same-site mean
%                            (z-space) — use for ROI/site comparisons
%      .t, .p, .df, .dz      paired-t result on spec_z and Cohen's dz
%      .mean_spec_z, .ci_spec_z
%      .metric, .tail
%
% :Examples:
% ::
%   % same > different within an ROI:
%   out = subject_level_crosssubject_effect(RSA_dp_corr);
%
%   % dpIns vs S1 cross-subject consistency, subject-level (paired across subjects):
%   o_dp = subject_level_crosssubject_effect(RSA_dp_corr);
%   o_s1 = subject_level_crosssubject_effect(RSA_s1_corr);
%   sub_dp = mean(o_dp.sitewise_same_z, 2, 'omitnan');   % one value/subject
%   sub_s1 = mean(o_s1.sitewise_same_z, 2, 'omitnan');
%   [~, p, ~, st] = ttest(sub_dp, sub_s1);   % df = nSub-1, NOT inflated
%
% :See also: permutation_test_site_specificity, get_crosssubject_site_effect
%
% ..
%    2026 Michael Sun. GPL v3.
% ..

p = inputParser;
p.addParameter('tail', 'right', @(x) any(strcmpi(x, {'right','both','left'})));
p.addParameter('doverbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
tail      = lower(p.Results.tail);
doverbose = logical(p.Results.doverbose);

nSub   = RSA.nSub;
nSites = RSA.nSites;
M      = RSA.M;
iscorr = strcmp(RSA.metric, 'correlation');

% Transform once: work in z-space for correlation, raw otherwise.
if iscorr
    Mt = rsa_fisher_z(M);
else
    Mt = M;
end

same_z = nan(nSub,1);
diff_z = nan(nSub,1);
sitewise_same_z = nan(nSub, nSites);

for a = 1:nSub
    is_a    = RSA.subject_idx == a;          % nodes belonging to anchor subject
    colsel  = RSA.subject_idx ~= a;          % partner nodes from other subjects

    same_vals = [];
    diff_vals = [];

    for b = 1:nSites
        % anchor node for site b
        a_node = is_a & (RSA.site_idx == b);
        if ~any(a_node), continue; end

        % same-site partners in other subjects (site b)
        same_partners = colsel & (RSA.site_idx == b);
        sv = Mt(a_node, same_partners);
        sv = sv(isfinite(sv));
        same_vals = [same_vals; sv(:)]; %#ok<AGROW>
        sitewise_same_z(a,b) = mean(sv, 'omitnan');

        % different-site partners in other subjects (site ~= b)
        diff_partners = colsel & (RSA.site_idx ~= b);
        dv = Mt(a_node, diff_partners);
        dv = dv(isfinite(dv));
        diff_vals = [diff_vals; dv(:)]; %#ok<AGROW>
    end

    same_z(a) = mean(same_vals, 'omitnan');
    diff_z(a) = mean(diff_vals, 'omitnan');
end

spec_z = same_z - diff_z;

% Paired test = one-sample t on the per-subject difference.
[~, pval, ci, st] = ttest(spec_z, 0, 'tail', tail);
dz = mean(spec_z, 'omitnan') / std(spec_z, 'omitnan');   % Cohen's dz (paired)

out = struct();
out.same_z          = same_z;
out.diff_z          = diff_z;
out.spec_z          = spec_z;
out.sitewise_same_z = sitewise_same_z;
out.t               = st.tstat;
out.p               = pval;
out.df              = st.df;
out.dz              = dz;
out.mean_spec_z     = mean(spec_z, 'omitnan');
out.ci_spec_z       = ci;
out.metric          = RSA.metric;
out.tail            = tail;

if iscorr
    out.same_r = tanh(same_z);
    out.diff_r = tanh(diff_z);
else
    out.same_r = same_z;   % already raw
    out.diff_r = diff_z;
end

if doverbose
    unit = 'z'; if ~iscorr, unit = RSA.metric; end
    fprintf('\nSubject-level cross-subject site specificity (metric=%s)\n', RSA.metric);
    fprintf('  same-site mean (r): %.4f | different-site mean (r): %.4f\n', ...
        mean(out.same_r,'omitnan'), mean(out.diff_r,'omitnan'));
    fprintf('  per-subject same-diff (%s): mean=%.4f, 95%% CI [%.4f %.4f]\n', ...
        unit, out.mean_spec_z, ci(1), ci(2));
    fprintf('  paired t(%d) = %.2f, p = %.4g (%s), Cohen''s dz = %.2f, N=%d subjects\n', ...
        out.df, out.t, out.p, tail, out.dz, nSub);
end

end
