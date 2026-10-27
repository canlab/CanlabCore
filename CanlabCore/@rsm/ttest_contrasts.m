function T = ttest_contrasts(obj, spec, varargin)
% ttest_contrasts  Run a battery of paired/one-sample t-tests across replicates.
%
% For each contrast in the spec, computes a per-replicate scalar via
% R.contrast() and runs ttest across replicates. Applies multiple-comparison
% correction (FDR / Bonferroni / none) across the battery.
%
% Usage
% -----
%   spec = { ...
%       % {name,         cells_A,         cells_B}
%       'within_hot',    'hot',           [];        % one-sample (cells_B empty)
%       'within_warm',   'warm',          [];
%       'hot_vs_warm',   'hot',           'warm';    % paired (within-hot vs within-warm)
%       'HI_vs_HW',      {'hot','imag'},  {'hot','warm'};
%       };
%   T = R.ttest_contrasts(spec, 'tail','right', 'correction','fdr');
%
% Returns a table with one row per contrast and columns:
%   Contrast, Mean_Diff, SE, t, df, p, FDR_P or Bonf_P, sig, Cohens_d
%
% Inputs
% ------
%   spec  N x 3 cell array: {name, cells_A, cells_B}
%         OR  N x 2 cell array (one-sample tests, cells_B implied empty)
%         Each cells_A / cells_B is a cells spec: name / {a,b} / index vector / [].
%
%   varargin:
%     'tail'        'both' (default) | 'right' | 'left'
%     'correction'  'fdr' (default) | 'bonferroni' | 'none'
%     'q'           threshold for sig flag (default 0.05)
%     'transform'   'auto' (default) | 'fisherz' | 'none'
%     'reduction'   'mean' (default) | 'median' | 'sum'

if builtin('numel', obj) > 1
    error('rsm:ttest_contrasts:nonScalar', 'expects a scalar rsm.');
end

p = inputParser;
p.addParameter('tail',       'both',  @(x) ischar(x) || isstring(x));
p.addParameter('correction', 'fdr',   @(x) ischar(x) || isstring(x));
p.addParameter('q',          0.05,    @isnumeric);
p.addParameter('transform',  'auto',  @(x) ischar(x) || isstring(x));
p.addParameter('reduction',  'mean',  @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opt = p.Results;

% Normalize spec
spec = normalize_spec(spec);
n = size(spec, 1);

names  = strings(n, 1);
diffs  = nan(n, 1);
ses    = nan(n, 1);
ts     = nan(n, 1);
dfs    = nan(n, 1);
ps     = nan(n, 1);
ds     = nan(n, 1);

for i = 1:n
    names(i) = string(spec{i, 1});
    A = spec{i, 2}; B = spec{i, 3};

    if isempty(B)
        % One-sample: test cells_A against zero
        v = obj.contrast(A, 'transform', opt.transform, 'reduction', opt.reduction);
        [~, p_i, ~, st] = ttest(v, 0, 'Tail', char(opt.tail));
        diffs(i) = mean(v, 'omitnan');
        ses(i) = std(v, 0, 'omitnan') / sqrt(sum(~isnan(v)));
        ds(i) = mean(v, 'omitnan') / std(v, 0, 'omitnan');
    else
        % Two-sample paired: test (cells_A - cells_B) against zero
        vA = obj.contrast(A, 'transform', opt.transform, 'reduction', opt.reduction);
        vB = obj.contrast(B, 'transform', opt.transform, 'reduction', opt.reduction);
        d = vA - vB;
        [~, p_i, ~, st] = ttest(d, 0, 'Tail', char(opt.tail));
        diffs(i) = mean(d, 'omitnan');
        ses(i) = std(d, 0, 'omitnan') / sqrt(sum(~isnan(d)));
        ds(i) = mean(d, 'omitnan') / std(d, 0, 'omitnan');
    end

    ts(i) = st.tstat;
    dfs(i) = st.df;
    ps(i) = p_i;
end

% Multiple comparison correction
correction = lower(char(opt.correction));
switch correction
    case 'fdr'
        if exist('FDR', 'file') == 2
            try
                p_thresh = FDR(ps, opt.q);
                if isempty(p_thresh) || isnan(p_thresh), p_thresh = 0; end
                corrected_p = ps;
                sig = ps <= p_thresh;
            catch
                % BH fallback
                [corrected_p, sig] = local_fdr_bh(ps, opt.q);
            end
        else
            [corrected_p, sig] = local_fdr_bh(ps, opt.q);
        end
        corrected_label = 'FDR_P';
    case 'bonferroni'
        corrected_p = min(ps * n, 1);
        sig = corrected_p < opt.q;
        corrected_label = 'Bonf_P';
    case 'none'
        corrected_p = ps;
        sig = ps < opt.q;
        % Avoid a duplicate 'P' column name -- use a distinct label.
        corrected_label = 'P_uncorr';
    otherwise
        error('rsm:ttest_contrasts:badCorrection', ...
            'correction must be ''fdr'', ''bonferroni'', or ''none''.');
end

% Build output table
T = table(names, diffs, ses, ts, dfs, ps, corrected_p, sig, ds, ...
    'VariableNames', {'Contrast', 'Mean_Diff', 'SE', 't', 'df', ...
                      'P', corrected_label, 'sig', 'Cohens_d'});

end


function spec = normalize_spec(spec)
% Accept either Nx2 (one-sample) or Nx3 (mixed) spec. Coerce to Nx3.
if size(spec, 2) == 2
    spec = [spec, cell(size(spec, 1), 1)];
end
if size(spec, 2) ~= 3
    error('rsm:ttest_contrasts:badSpec', ...
        'spec must be Nx2 (one-sample) or Nx3 ({name, A, B}).');
end
end


function [corrected_p, sig] = local_fdr_bh(ps, q)
% Benjamini-Hochberg stock fallback
n = numel(ps);
[sorted_p, sort_idx] = sort(ps);
thresh = (1:n)' / n * q;
below = sorted_p(:) <= thresh;
if any(below)
    k = find(below, 1, 'last');
    p_thresh = sorted_p(k);
else
    p_thresh = 0;
end
% Storey-style adjusted p; simpler version: min(p * n / rank, 1)
ranks = zeros(n, 1); ranks(sort_idx) = 1:n;
corrected_p = min(ps(:) .* n ./ ranks, 1);
sig = ps(:) <= p_thresh;
end
