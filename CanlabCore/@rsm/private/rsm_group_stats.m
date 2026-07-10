function G = rsm_group_stats(V, correction, nperm)
% rsm_group_stats  One-sample group stats over subjects for RSA count methods.
%
% Thin wrapper used by @rsm/count_map, count_models, count_regions to get a
% group effect + p-value for each cell/model/region across subjects. When the
% shared permutation engine hrf_group_stats (HRF_Toolbox_Pipeline) is on the
% path it is used, so RSA inference matches the rest of the CANlab HRF API;
% otherwise a built-in one-sample t (with FDR/Bonferroni) is used, keeping the
% RSA toolbox self-contained.
%
% :Inputs:
%   **V:** [nCells x N] (rows = cells/models/regions, cols = subjects) or a
%          [1 x N] vector for a single cell.
%   **correction:** 'fdr' (default) | 'permutation' | 'bonferroni' | 'none'.
%   **nperm:** permutations when hrf_group_stats does the work. Default 5000.
%
% :Output:
%   **G:** struct with .est, .t, .p, .p_corr, .sig (one row per cell).

if nargin < 2 || isempty(correction), correction = 'fdr'; end
if nargin < 3 || isempty(nperm), nperm = 5000; end
V = double(V);
if isvector(V), V = V(:)'; end

if exist('hrf_group_stats', 'file') == 2
    G = hrf_group_stats(V, 'Correction', correction, 'Nperm', nperm);
    return
end

% ---- built-in fallback: one-sample t across subjects (cols) -------------
[nc, N] = size(V);
est = mean(V, 2, 'omitnan');
se  = std(V, 0, 2, 'omitnan') / sqrt(max(N, 1));
t   = est ./ se; t(se == 0) = 0;
p   = 2 * (1 - local_tcdf(abs(t), N - 1));
switch lower(char(correction))
    case 'none',       pc = p;
    case 'bonferroni', pc = min(p * nc, 1);
    otherwise,         pc = local_fdr_bh(p);      % 'fdr' (and 'permutation' w/o engine)
end
G = struct('est', est, 't', t, 'p', p, 'p_corr', pc, 'sig', pc < 0.05);
end


function pcdf = local_tcdf(tval, df)
if df <= 0, pcdf = 0.5 * ones(size(tval)); return; end
x = df ./ (df + tval .^ 2);
pcdf = 1 - 0.5 * betainc(x, df / 2, 0.5);
end


function q = local_fdr_bh(p)
% Benjamini-Hochberg adjusted p-values.
p = p(:); m = numel(p);
[ps, ord] = sort(p);
q = ps .* m ./ (1:m)';
q = min(1, flipud(cummin(flipud(q))));
out = nan(m, 1); out(ord) = q; q = out;
end
