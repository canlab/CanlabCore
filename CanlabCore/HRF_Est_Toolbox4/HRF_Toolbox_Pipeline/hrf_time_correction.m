function C = hrf_time_correction(D, varargin)
%HRF_TIME_CORRECTION Multiple-comparison correction across a peristimulus timecourse.
%
% Thin adapter over hrf_group_stats for the [nSubj x nTime] per-subject
% contrast matrices used by the HRF timecourse tests (hrf_time_unfolding_stats,
% hrf_2x2_study_score_stats, hrf_compare_conditions). It gives the per-lag
% inference a consistent, small-n-appropriate correction -- sign-flip / label
% max-|t| FWER, temporal cluster-mass, or FDR -- instead of the uncorrected
% per-timepoint t-tests those functions computed on their own. This is the
% correction the n~7-11 per-lag HRF summaries need.
%
% :Usage:
% ::
%     C = hrf_time_correction(D, 'Correction','cluster')                % one-sample vs 0
%     C = hrf_time_correction(D, 'Group', g, 'Correction','permutation')% two-sample diff
%
% :Inputs:
%   **D:** [nSubj x nTime] per-subject timecourse. One-sample tests each
%          timepoint's mean against 0 (D is usually a within-subject contrast).
%
% :Optional Inputs:
%   **'Correction':** 'none' (default) | 'fdr' | 'permutation' | 'cluster'.
%   **'Group':** [nSubj] two-group labels; switches to a two-sample test of the
%                group difference at each timepoint.
%   **'Nperm':** permutations when the exact set is too large. Default 5000.
%   **'Alpha':** corrected-p / FWER threshold. Default 0.05.
%
% :Output:
%   **C:** struct with column vectors [nTime x 1] .t, .p, .p_corr, .sig, plus
%          .correction, .design, .n.
%
% See also: hrf_group_stats, hrf_time_unfolding_stats, hrf_2x2_study_score_stats,
%           hrf_compare_conditions.

p = inputParser;
p.addRequired('D', @isnumeric);
p.addParameter('Correction', 'none', @(x) ischar(x) || isstring(x));
p.addParameter('Group', [], @(x) isempty(x) || isvector(x) || iscell(x));
p.addParameter('Nperm', 5000, @(x) isscalar(x) && x >= 100);
p.addParameter('Alpha', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
p.parse(D, varargin{:});
o = p.Results;

[nsub, ntime] = size(D);
corr = lower(char(o.Correction));

% [nTime x 1 x nSubj] so hrf_group_stats treats time as the tested family and
% the last dim as subjects; the singleton keeps it 2-D for cluster contiguity.
data = reshape(D.', [ntime, 1, nsub]);

args = {'Correction', corr, 'Nperm', o.Nperm, 'Threshold', o.Alpha};
if ~isempty(o.Group)
    args = [args, {'Design', 'twosample', 'Group', o.Group}];
end
if any(strcmp(corr, {'cluster', 'tcluster'}))
    args = [args, {'ClusterDim', 1}];
end

G = hrf_group_stats(data, args{:});

C = struct();
C.t = G.t(:);
C.p = G.p(:);
C.p_corr = G.p_corr(:);
C.sig = G.sig(:);
C.correction = corr;
C.design = G.design;
C.n = nsub;
end
