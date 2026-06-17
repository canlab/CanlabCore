function fig = hrf_misspec_report(M, varargin)
% One-figure dashboard for the HRF misspecification pipeline + auto next-steps.
%
% Reduces the curve-misspec table (hrf_misspec_metrics) and, optionally, the
% residual table (hrf_residual_diagnostics) to a four-panel figure you can
% eyeball in seconds:
%   A  misspec_r2 by region/source   - where canonical fits (green) vs fails (red)
%   B  peak_lag_bias by region        - which HRFs are delayed (+) / advanced (-)
%   C  residual model comparison      - which model leaves task-locked structure
%   D  DIAGNOSTICS & NEXT STEPS       - auto-flags (noisy curves, residual
%                                       autocorrelation, uncomputed power-loss, ...)
%
% :Usage:
% ::
%     hrf_misspec_report(M)                          % curve half only
%     hrf_misspec_report(M, 'Residuals', R)          % + residual half
%     hrf_misspec_report(M, 'Condition','rest_stim', 'TopN',12, 'Save','report.png')
%
% :Inputs:
%   **M:** table from hrf_misspec_metrics (per-subject OR GroupCurveFirst).
%
% :Optional Inputs:
%   **'Residuals':** table from hrf_residual_diagnostics. Default none.
%   **'Condition':** restrict M before reducing. char/string/cellstr; glob
%                    wildcards ok ('*heat*', 'rest_stim_ttl_?'); multiple
%                    patterns match any (and pool into the per-region mean).
%                    Default '' = all conditions.
%   **'TopN':**      how many best/worst regions to bar. Default 14.
%   **'GoodR2':**    misspec_r2 above this counts as well-described. Default 0.5.
%   **'Title':**     figure title. Default 'HRF misspecification report'.
%   **'Save':**      path to save a PNG (also returns the handle). Default ''.
%
% :Output:  **fig** - the figure handle.
%
% See also: hrf_misspec_metrics, hrf_residual_diagnostics, hrf_curve_summary_groupstats.

p = inputParser;
p.addRequired('M', @istable);
p.addParameter('Residuals', table(), @istable);
p.addParameter('Condition', '', @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('TopN', 14, @(x) isscalar(x) && x >= 1);
p.addParameter('GoodR2', 0.5, @(x) isscalar(x));
p.addParameter('Title', 'HRF misspecification report', @(x) ischar(x) || isstring(x));
p.addParameter('Save', '', @(x) ischar(x) || isstring(x));
p.parse(M, varargin{:});
opts = p.Results;
R = opts.Residuals;

% ---- reduce curve table to per-region means -----------------------------
Mc = M;
pats = string(opts.Condition);
pats = pats(strlength(strtrim(pats)) > 0);
if ~isempty(pats) && any(strcmp('condition', M.Properties.VariableNames))
    Mc = M(local_match_conditions(M.condition, pats), :);   % glob wildcards ok
end
if isempty(Mc) || height(Mc) == 0
    avail = '';
    if any(strcmp('condition', M.Properties.VariableNames))
        avail = strjoin(unique(cellstr(string(M.condition)), 'stable'), ', ');
    end
    error('hrf_misspec_report:NoRows', ...
        'No rows in M after Condition filter {%s}. Available conditions: %s', ...
        strjoin(cellstr(pats), ', '), avail);
end
reg = string(Mc.source_name);
[g, ureg] = findgroups(reg);
r2  = splitapply(@(x) mean(x, 'omitnan'), Mc.misspec_r2, g);
plb = splitapply(@(x) mean(x, 'omitnan'), Mc.peak_lag_bias_seconds, g);
nsub = splitapply(@numel, Mc.misspec_r2, g);

fig = figure('Color', 'w', 'Position', [80 80 1280 860], 'Name', char(opts.Title));
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, char(opts.Title), 'FontWeight', 'bold', 'Interpreter', 'none');

% ---- Panel A: misspec_r2 ------------------------------------------------
nexttile(tl);
local_ranked_barh(ureg, r2, opts.TopN, [0 opts.GoodR2], ...
    'misspec_r2  (1 = canonical fits; \leq0 = canonical fails)', true);

% ---- Panel B: peak_lag_bias --------------------------------------------
nexttile(tl);
local_ranked_barh(ureg, plb, opts.TopN, 0, ...
    'peak\_lag\_bias (s)   (+ delayed vs canonical, - advanced)', false);

% ---- Panel C: residual model comparison --------------------------------
ax = nexttile(tl);
[resid_summary, have_resid] = local_resid_summary(R);
if have_resid
    local_resid_panel(ax, resid_summary);
else
    axis(ax, 'off');
    text(ax, 0.5, 0.5, {'(no residual table supplied)', ...
        'pass ''Residuals'', hrf\_residual\_diagnostics(study)'}, ...
        'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);
end

% ---- Panel D: diagnostics & next steps ---------------------------------
ax = nexttile(tl);
axis(ax, 'off');
lines = local_diagnostics(r2, plb, nsub, opts.GoodR2, resid_summary, have_resid, Mc);
text(ax, 0.0, 1.0, lines, 'VerticalAlignment', 'top', 'FontName', 'FixedWidth', ...
    'FontSize', 9, 'Interpreter', 'none');
title(ax, 'Diagnostics & next steps', 'FontWeight', 'bold');

if ~isempty(char(opts.Save))
    try
        exportgraphics(fig, char(opts.Save), 'Resolution', 130);
    catch
        saveas(fig, char(opts.Save));
    end
end
end


% =========================================================================
function local_ranked_barh(labels, vals, topN, ref_lines, ttl, color_by_value)
ok = isfinite(vals);
labels = labels(ok); vals = vals(ok);
[vals, ix] = sort(vals, 'ascend'); labels = labels(ix);
n = numel(vals);
if n > 2 * topN
    keep = [1:topN, (n - topN + 1):n];
    vals = vals(keep); labels = labels(keep);
end
b = barh(vals, 'FaceColor', 'flat');
if color_by_value
    b.CData = local_value_colors(vals);
else
    c = repmat([0.30 0.45 0.80], numel(vals), 1);   % advanced = blue
    c(vals > 0, :) = repmat([0.80 0.30 0.30], sum(vals > 0), 1); % delayed = red
    b.CData = c;
end
set(gca, 'YTick', 1:numel(vals), 'YTickLabel', labels, 'TickLabelInterpreter', 'none', 'FontSize', 7.5);
ylim([0.5 numel(vals) + 0.5]);
grid on; box off;
xlabel(ttl, 'Interpreter', 'tex', 'FontSize', 9);
for r = ref_lines(:)'
    xline(r, '--', 'Color', [0.5 0.5 0.5]);
end
end

function C = local_value_colors(v)
% red (low/negative) -> yellow -> green (high), clamped to [0,1] of value.
t = max(min(v, 1), 0);
C = [1 - t, 0.55 + 0.30 * t, 0.25 * ones(size(t))];
C(v < 0, :) = repmat([0.70 0.10 0.10], sum(v < 0), 1);   % negative = deep red
end

function [S, ok] = local_resid_summary(R)
S = struct('model', {{}}, 'mis_modeling_p', [], 'durbin_watson', [], 'power_loss', []);
ok = false;
if isempty(R) || height(R) == 0 || ~any(strcmp('model', R.Properties.VariableNames)), return; end
mdl = string(R.model);
[g, um] = findgroups(mdl);
S.model = cellstr(um);
S.mis_modeling_p = splitapply(@(x) median(x, 'omitnan'), R.mis_modeling_p, g);
S.durbin_watson  = splitapply(@(x) median(x, 'omitnan'), R.durbin_watson, g);
S.power_loss     = splitapply(@(x) median(x, 'omitnan'), R.power_loss, g);
ok = true;
end

function local_resid_panel(ax, S)
nm = S.model;
dw = S.durbin_watson;
b = bar(ax, dw, 'FaceColor', [0.45 0.55 0.70]);  %#ok<NASGU>
set(ax, 'XTick', 1:numel(nm), 'XTickLabel', nm, 'TickLabelInterpreter', 'none');
yline(ax, 2, '--', 'white (DW=2)', 'Color', [0.2 0.2 0.2]);
ylabel(ax, 'Durbin-Watson (residual)'); grid(ax, 'on'); box(ax, 'off');
% annotate mis_modeling_p / power_loss per model under each bar
for i = 1:numel(nm)
    text(ax, i, max(dw) * 0.06 + 0.02, ...
        sprintf('p=%.1g\nPL=%.2g', S.mis_modeling_p(i), S.power_loss(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 7.5, 'Color', [0.25 0.25 0.25]);
end
title(ax, 'Residual structure by model  (lower DW / smaller p = worse)', 'FontSize', 9);
end

function lines = local_diagnostics(r2, plb, nsub, goodR2, S, have_resid, Mc)
L = {};
maxr2 = max(r2);
frac_neg = mean(r2 < 0);
frac_good = mean(r2 >= goodR2);
n_curves_per_region = round(median(nsub));
is_group = any(strcmpi(string(Mc.Properties.VariableNames), 'subject')) && ...
    all(string(Mc.subject) == "group");

L{end+1} = sprintf('CURVE FIT (n regions = %d, %s)', numel(r2), ...
    local_tern(is_group, 'GroupCurveFirst', sprintf('~%d curves/region', n_curves_per_region)));
L{end+1} = sprintf('  best misspec_r2 = %+.2f | %.0f%% well-fit (>=%.2f) | %.0f%% canonical-fails (<0)', ...
    maxr2, 100 * frac_good, goodR2, 100 * frac_neg);
if maxr2 < 0.4
    L{end+1} = '  [!] best fit is poor -> curves look NOISE-DOMINATED.';
    if ~is_group
        L{end+1} = '      try GroupCurveFirst=true; and/or collapse sparse';
        L{end+1} = '      conditions (avoid *_ttl_* splits) when fitting.';
    else
        L{end+1} = '      group curve already used -> the FIT is noisy:';
        L{end+1} = '      collapse sparse conditions / model the block.';
    end
else
    L{end+1} = '  [ok] canonical describes a subset of regions well.';
end
[mxbias, mxi] = max(plb);
L{end+1} = sprintf('  largest +lag bias: %+.1f s  (most delayed HRF)', mxbias); %#ok<NASGU>
L{end+1} = '';
L{end+1} = 'RESIDUALS';
if have_resid
    [~, worst] = min(S.durbin_watson);   % furthest from white
    L{end+1} = sprintf('  median DW: %s', strjoin(arrayfun(@(i) sprintf('%s=%.2f', S.model{i}, S.durbin_watson(i)), 1:numel(S.model), 'uni', 0), '  '));
    if min(S.durbin_watson) < 1.5
        L{end+1} = '  [!] DW << 2 -> strong residual AUTOCORRELATION.';
        L{end+1} = '      mis_modeling_p then reflects autocorrelation, not';
        L{end+1} = '      HRF fit -> WHITEN the 1D fits before trusting it.';
    end
    if all(S.power_loss == 0 | isnan(S.power_loss))
        L{end+1} = '  [!] power_loss all 0/NaN -> not computed at fit time.';
        L{end+1} = '      pass Recompute=true (needs the fit inputs).';
    end
    L{end+1} = sprintf('  most residual structure: model "%s"', S.model{worst});
else
    L{end+1} = '  (supply Residuals to compare models)';
end
L{end+1} = '';
L{end+1} = 'NEXT STEPS';
if maxr2 < 0.4
    L{end+1} = '  1) fix the condition (collapse/block), re-fit, re-score';
    L{end+1} = '  2) re-run this report; expect best misspec_r2 -> 0.7-0.9';
else
    L{end+1} = '  1) inspect the red (low-r2) regions'' HRF curves vs canonical';
    L{end+1} = '  2) report peak_lag_bias for the significant delayed regions';
end
lines = L(:);
end

function s = local_tern(c, a, b)
if c, s = a; else, s = b; end
end

function mask = local_match_conditions(conds, patterns)
% Glob-match condition strings against one or more wildcard patterns
% (e.g. '*heat*', 'rest_stim_ttl_?'). Exact strings still match exactly.
conds = cellstr(string(conds));
patterns = cellstr(string(patterns));
mask = false(numel(conds), 1);
for p = 1:numel(patterns)
    pat = strtrim(patterns{p});
    if isempty(pat), continue; end
    if any(pat == '*' | pat == '?')
        rx = ['^', regexptranslate('wildcard', pat), '$'];
        hit = ~cellfun('isempty', regexp(conds, rx, 'once'));
    else
        hit = strcmp(conds, pat);
    end
    mask = mask | hit(:);
end
end
