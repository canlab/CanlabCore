function out = hrf_animate_wordcloud(source, varargin)
%HRF_ANIMATE_WORDCLOUD Animated term wordcloud of HRF map-scores over lags.
%
% Visualizes how term/map associations EVOLVE across peristimulus time
% (HRF lags): one wordcloud frame per lag, each term sized by its score at
% that lag and coloured by sign. Built for Neurosynth term-map scores
% (map_neurosynth_<term> columns) but works for any map_<set>_* score block.
%
% Terms are placed at FIXED positions (a spiral ordered by overall
% importance) so only the font size/colour animate -- you can actually see a
% term grow and fade as the haemodynamic response unfolds, instead of the
% layout jumping every frame (as MATLAB's wordcloud would).
%
% :Usage:
% ::
%     hrf_animate_wordcloud(score_csv, 'Condition','heat', 'OutputFile','terms.mp4')
%     hrf_animate_wordcloud(input_table, 'Set','neurosynth', 'TopN',40)  % group mean
%     hrf_animate_wordcloud({lf,obs}, 'Set','neurosynth', 'Model','sfir', ...
%         'Condition','*heat*', 'OutputFile','pooled_terms.mp4')   % POOL dirs
%     out = hrf_animate_wordcloud(struct('scores',M,'terms',t,'lags',L))   % direct
%
% :Inputs:
%   **source:** any of -- a score CSV path; an input_table
%             (subject/*_scores_file rows); an output DIRECTORY or a CELL of
%             directories/input tables (collected and POOLED into one group
%             mean across every subject of every dir -- same subject id across
%             dirs combines); or a struct with fields .scores [nLag x nTerm],
%             .terms (1 x nTerm), .lags (1 x nLag).
%
% :Optional Inputs:
%   **'Set':**        map set token (default 'neurosynth'); selects
%                     map_<set>_<term> columns.
%   **'Condition':**  condition (glob ok) whose curve to animate. Default '' =
%                     first condition present (errors if several and ambiguous).
%   **'Model'/'Object':** for input_table iteration. Default 'sfir'/'beta'.
%   **'Threshold':**  significance cutoff for selecting/displaying terms, via a
%                     one-sample t across subjects at each (term,lag). Default
%                     0.05. A term is SHOWN at the lags where its association is
%                     significant (it lights up over the HRF), and a term is
%                     selected only if significant at >=1 lag. Needs a group
%                     (>=2 subjects); with a single CSV/struct there is no group
%                     so it falls back to top-N by |score|.
%   **'Correction':** 'fdr' (default, BH across terms within each lag) or
%                     'none'.
%   **'Persist':**    false (default) -- a selected term is drawn only at lags
%                     where it is significant. true -- always draw selected
%                     terms, greyed at sub-threshold lags.
%   **'TopN':**       LEGIBILITY CAP on how many terms are shown (default 60),
%                     applied AFTER the statistical selection, keeping the most
%                     significant. (Only when there is no group does it act as
%                     a pure top-N by |score|.)
%   **'SizeBy':**     'abs' (default) | 'pos' | 'signed' -- which magnitude
%                     drives font size (abs = both directions grow).
%   **'FrameRate':**  fps of the movie. Default 4.
%   **'OutputFile':** '.mp4'/'.avi'/'.gif' to write; '' = just show. Default ''.
%   **'FontRange':**  [min max] point sizes. Default [9 46].
%   **'Title':**      title prefix. Default 'Neurosynth terms over the HRF'.
%   **'Verbose'/'doverbose':** chatter. Default true.
%
% :Outputs:
%   **out:** struct with .scores/.t/.p/.sig [nLag x nTerm_kept], .terms, .lags,
%            .pos (fixed [k x 2] layout), .nsubj, .selection (label), .file.
%
% :Examples:
% ::
%     % after rescuing neurosynth scores (see hrf_score_one_prefix append mode):
%     IT = hrf_collect_wholebrain_outputs(out_dir);
%     hrf_animate_wordcloud(IT, 'Set','neurosynth', 'Condition','*heat*', ...
%         'Model','sfir', 'TopN',40, 'OutputFile', fullfile(out_dir,'heat_terms.mp4'));
%
% See also: hrf_score_one_prefix, hrf_apply_maps_to_wholebrain,
%           neurosynth_lexical_plot, hrf_misspec_metrics.

p = inputParser;
p.addRequired('source');
p.addParameter('Set', 'neurosynth', @(x) ischar(x) || isstring(x));
p.addParameter('Condition', '', @(x) ischar(x) || isstring(x) || iscell(x));
p.addParameter('Model', 'sfir', @(x) ischar(x) || isstring(x));
p.addParameter('Object', 'beta', @(x) ischar(x) || isstring(x));
p.addParameter('Threshold', 0.05, @(x) isscalar(x) && x > 0 && x <= 1);
p.addParameter('Correction', 'fdr', @(x) ischar(x) || isstring(x));
p.addParameter('Persist', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('TopN', 60, @(x) isscalar(x) && x >= 1);
p.addParameter('SizeBy', 'abs', @(x) ischar(x) || isstring(x));
p.addParameter('FrameRate', 4, @(x) isscalar(x) && x > 0);
p.addParameter('OutputFile', '', @(x) ischar(x) || isstring(x));
p.addParameter('FontRange', [9 46], @(x) isnumeric(x) && numel(x) == 2);
p.addParameter('Title', 'Neurosynth terms over the HRF', @(x) ischar(x) || isstring(x));
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('doverbose', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
p.parse(source, varargin{:});
opts = p.Results;
verbose = logical(opts.Verbose);
if ~isempty(opts.doverbose), verbose = logical(opts.doverbose); end

S = local_get_term_stats(source, opts);          % .scores/.t/.p [nLag x nTerm], .nsubj
if isempty(S.scores), error('hrf_animate_wordcloud:NoData', 'No map_%s_* scores found for the requested condition.', char(opts.Set)); end
scores = S.scores; terms = S.terms; lags = S.lags; condlabel = S.condlabel;
[nLag, nT] = size(scores);
mag = local_size_metric(scores, opts.SizeBy);

% ---- significance: which (term,lag) associations are non-zero across subjects
have_stats = isfield(S, 'p') && S.nsubj >= 2 && any(isfinite(S.p(:)));
thr = opts.Threshold; corr = lower(char(opts.Correction));
if have_stats
    sig = local_sig_mask(S.p, thr, corr);        % [nLag x nT] logical
    sel_note = sprintf('p<%.3g %s, group n=%d', thr, upper(corr), S.nsubj);
else
    sig = true(nLag, nT);
    if verbose
        warning('hrf_animate_wordcloud:NoGroupStats', ...
            'No across-subject statistics (need >=2 subjects); selecting top %d terms by |score|.', opts.TopN);
    end
    sel_note = sprintf('top %d by |score| (no group stats)', opts.TopN);
end

% ---- select terms (significant at any lag), rank, cap at TopN -----------
passed = any(sig, 1);
gate = have_stats;                                % per-frame significance gating
if ~any(passed)
    warning('hrf_animate_wordcloud:NoneSignificant', ...
        'No term significant at p<%.3g (%s); showing top %d by |score| instead.', thr, upper(corr), opts.TopN);
    passed = true(1, nT); sig = true(nLag, nT); gate = false;
    sel_note = sprintf('top %d by |score| (none reached p<%.3g)', opts.TopN, thr);
end
if have_stats && any(any(sig))
    rankval = min(S.p, [], 1, 'omitnan'); rankval(~passed) = Inf;
    [~, ord] = sort(rankval, 'ascend');
else
    imp = max(mag, [], 1, 'omitnan'); imp(~passed) = -Inf;
    [~, ord] = sort(imp, 'descend');
end
keepN = min(opts.TopN, sum(passed));
keep = ord(1:keepN);
terms = terms(keep); scores = scores(:, keep); mag = mag(:, keep); sig = sig(:, keep);

gmax = max(mag(:)); if gmax == 0 || ~isfinite(gmax), gmax = 1; end
clim = max(abs(scores(:))); if clim == 0 || ~isfinite(clim), clim = 1; end
fr = opts.FontRange;

fig = figure('Color', 'w', 'Position', [80 80 1000 760], 'Name', char(opts.Title));
ax = axes(fig, 'Position', [0.02 0.04 0.96 0.90]); axis(ax, [0 1 0 1]); axis(ax, 'off'); hold(ax, 'on');

% Non-overlapping wordcloud packing computed ONCE at each word's MAXIMUM font
% size (its peak over lags). Every later frame only shrinks a word, so the
% packing stays overlap-free and positions are stable across the movie -- the
% tight look of a static wordcloud, but animated.
peakmag = max(mag, [], 1, 'omitnan');
pos = local_wordcloud_layout(ax, terms, peakmag, gmax, fr);

vw = local_open_video(opts.OutputFile, opts.FrameRate);
persist = logical(opts.Persist);
for li = 1:nLag
    cla(ax); axis(ax, [0 1 0 1]); axis(ax, 'off');
    for ti = 1:numel(terms)
        s = scores(li, ti); m = mag(li, ti);
        if ~isfinite(m) || m <= 0, continue; end
        % gate by per-lag significance (a term lights up only where it is
        % significant), unless Persist or no group stats.
        if gate && ~persist && ~sig(li, ti), continue; end
        fs = fr(1) + (fr(2) - fr(1)) * (m / gmax);
        if gate && persist && ~sig(li, ti)
            col = [0.7 0.7 0.7];                  % shown but greyed (sub-threshold)
        else
            col = local_sign_color(s, clim);
        end
        text(ax, pos(ti, 1), pos(ti, 2), char(terms{ti}), 'FontSize', fs, ...
            'Color', col, 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'Interpreter', 'none', 'Clipping', 'on');
    end
    title(ax, sprintf('%s   —   %s   t = %.1f s   [%s]', char(opts.Title), condlabel, lags(li), sel_note), ...
        'Interpreter', 'none', 'FontSize', 11);
    drawnow;
    vw = local_write_frame(vw, fig);
end
vw = local_close_video(vw);

out = struct('scores', scores, 'terms', {terms}, 'lags', lags(:)', 'pos', pos, ...
    't', S.t(:, keep), 'p', S.p(:, keep), 'sig', sig, 'nsubj', S.nsubj, ...
    'selection', sel_note, 'file', local_video_path(vw, opts.OutputFile));
if verbose
    fprintf('hrf_animate_wordcloud: %d terms shown (%s) x %d lags, condition %s%s\n', ...
        numel(terms), sel_note, nLag, condlabel, local_file_note(out.file));
end
end


% =========================================================================
function S = local_get_term_stats(source, opts)
% Returns S.scores/.t/.p [nLag x nTerm], .terms, .lags, .condlabel, .nsubj.
% Group statistics (across-subject one-sample t) are produced for the
% input_table / dir / cell sources; single CSV and struct sources have no
% group (t/p NaN -> caller falls back to top-N by |score|).
prefix = ['map_', char(opts.Set), '_'];
S = struct('scores', [], 't', [], 'p', [], 'terms', {{}}, 'lags', [], 'condlabel', '', 'nsubj', 0);
if isstruct(source) && isfield(source, 'scores')
    S.scores = source.scores; S.terms = cellstr(string(source.terms));
    S.lags = source.lags(:)'; S.condlabel = 'curve';
    S.t = nan(size(S.scores)); S.p = nan(size(S.scores)); S.nsubj = 1;
    return
end
is_dir = (ischar(source) || isstring(source)) && isfolder(char(string(source)));
if iscell(source) || is_dir
    S = local_stats_from_input_table(local_pool_input_table(source), prefix, opts);
    return
end
if (ischar(source) || isstring(source)) && endsWith(string(source), '.csv')
    [M, terms, lags, cl] = local_matrix_from_csv(char(source), prefix, opts);
    S.scores = M; S.terms = terms; S.lags = lags; S.condlabel = cl;
    S.t = nan(size(M)); S.p = nan(size(M)); S.nsubj = 1;
    return
end
if istable(source)
    S = local_stats_from_input_table(source, prefix, opts);
    return
end
error('hrf_animate_wordcloud:Source', ...
    'source must be a score CSV, an input_table, a struct, an output dir, or a cell of dirs/tables.');
end


function IT = local_pool_input_table(source)
% Concatenate the collection tables of one or more output dirs (or pre-built
% input tables) on their common columns. Same subject id across dirs pools
% that subject's score files in the downstream group mean.
if iscell(source), items = source(:)'; else, items = {source}; end
IT = table();
for i = 1:numel(items)
    it = items{i};
    if istable(it)
        Ti = it;
    else
        Ti = hrf_collect_wholebrain_outputs(char(string(it)));
    end
    if isempty(Ti) || height(Ti) == 0, continue; end
    if isempty(IT) || height(IT) == 0
        IT = Ti;
    else
        c = intersect(IT.Properties.VariableNames, Ti.Properties.VariableNames, 'stable');
        IT = [IT(:, c); Ti(:, c)]; %#ok<AGROW>
    end
end
if isempty(IT) || height(IT) == 0
    error('hrf_animate_wordcloud:NoRecords', 'No score records collected from the given dir(s).');
end
end


function [M, terms, lags, condlabel] = local_matrix_from_csv(csv, prefix, opts)
T = readtable(csv, 'TextType', 'string');
[M, terms, lags, condlabel] = local_matrix_from_table(T, prefix, opts);
end


function S = local_stats_from_input_table(IT, prefix, opts)
% Per-SUBJECT term matrices (average that subject's runs/files), stacked, then
% group mean + one-sample t across subjects (the across-subject inference that
% the statistical threshold uses).
S = struct('scores', [], 't', [], 'p', [], 'terms', {{}}, 'lags', [], 'condlabel', '', 'nsubj', 0);
vars = IT.Properties.VariableNames;
file_col = '';
if strcmpi(char(opts.Object), 't') && any(strcmp('t_scores_file', vars))
    file_col = 't_scores_file';
elseif any(strcmp('beta_scores_file', vars))
    file_col = 'beta_scores_file';
end
if isempty(file_col), error('hrf_animate_wordcloud:NoScoreFiles', 'input_table lacks *_scores_file columns.'); end
model = char(opts.Model);
if any(strcmp('subject', vars)), subjects = cellstr(string(IT.subject));
else, subjects = arrayfun(@(i) sprintf('row%d', i), 1:height(IT), 'uni', 0); end
usub = unique(subjects, 'stable');

ref_terms = {}; ref_lags = []; condlabel = '';
parts = {}; used = {};
for s = 1:numel(usub)
    rows = find(strcmp(subjects, usub{s}));
    acc = []; nf = 0;
    for ii = rows(:)'
        if any(strcmp('model', vars)) && ~strcmpi(char(string(IT.model(ii))), model), continue; end
        f = char(string(IT.(file_col)(ii)));
        if isempty(f) || exist(f, 'file') ~= 2, continue; end
        [Mi, ti, li, cl] = local_matrix_from_table(readtable(f, 'TextType', 'string'), prefix, opts);
        if isempty(Mi), continue; end
        if isempty(ref_terms), ref_terms = ti; ref_lags = li; condlabel = cl; end
        if ~isequal(ti, ref_terms) || ~isequal(li, ref_lags)
            warning('hrf_animate_wordcloud:GridMismatch', 'Skipping a score file whose term/lag grid differs.');
            continue
        end
        if isempty(acc), acc = Mi; nf = 1; else, acc = acc + Mi; nf = nf + 1; end
    end
    if nf == 0, continue; end
    parts{end + 1} = acc / nf; used{end + 1} = usub{s}; %#ok<AGROW>
end
if isempty(parts), return; end

stack = cat(3, parts{:});
n = size(stack, 3);
m = mean(stack, 3, 'omitnan');
S.scores = m; S.terms = ref_terms; S.lags = ref_lags; S.condlabel = condlabel; S.nsubj = n;
if n >= 2
    se = std(stack, 0, 3, 'omitnan') / sqrt(n);
    t = m ./ se; t(se == 0) = 0;
    S.t = t; S.p = 2 * (1 - local_tcdf(abs(t), n - 1));
else
    S.t = nan(size(m)); S.p = nan(size(m));
end
end


function sig = local_sig_mask(p, thr, corr)
% Per-lag significance: at each lag, test terms; FDR-correct across terms
% within that lag (so each frame's "significant terms" is FDR-controlled).
[L, T] = size(p);
sig = false(L, T);
for li = 1:L
    pl = p(li, :);
    valid = isfinite(pl);
    if ~any(valid), continue; end
    if strcmp(corr, 'fdr')
        padj = local_fdr_row(pl(valid));
        row = false(1, T); row(valid) = padj <= thr; sig(li, :) = row;
    else
        sig(li, valid) = pl(valid) < thr;
    end
end
end

function padj = local_fdr_row(p)
% Benjamini-Hochberg adjusted p-values for a vector.
p = p(:)'; m = numel(p);
[sp, ord] = sort(p);
adj = min(1, fliplr(cummin(fliplr(sp .* m ./ (1:m)))));
padj = nan(1, m); padj(ord) = adj;
end

function pcdf = local_tcdf(tval, df)
% Student-t CDF via regularized incomplete beta (no Stats toolbox needed).
x = df ./ (df + tval .^ 2);
pcdf = 1 - 0.5 * betainc(x, df / 2, 0.5);
end


function [M, terms, lags, condlabel] = local_matrix_from_table(T, prefix, opts)
M = []; terms = {}; lags = []; condlabel = '';
v = T.Properties.VariableNames;
cols = v(startsWith(v, prefix) & ~endsWith(v, '_se'));
if isempty(cols), return; end
if any(strcmp('lag_seconds', v)), lagcol = 'lag_seconds'; else, lagcol = 'lag_index'; end
if ~any(strcmp('condition', v)) || ~any(strcmp(lagcol, v)), return; end
cond = string(T.condition);
cmask = local_cond_pick(cond, opts.Condition);
condlabel = char(local_cond_label(cond, cmask));
lg = double(T.(lagcol));
[ul, ~, gi] = unique(lg(cmask), 'stable');
[lags, so] = sort(ul);
terms = cellfun(@(c) local_term_name(c, prefix), cols, 'uni', 0);
M = zeros(numel(lags), numel(cols));
for j = 1:numel(cols)
    y = double(T.(cols{j})(cmask));
    m = accumarray(gi, y, [], @(x) mean(x, 'omitnan'));
    M(:, j) = m(so);
end
end


function name = local_term_name(col, prefix)
name = col(numel(prefix) + 1:end);
end


function mask = local_cond_pick(cond, want)
want = cellstr(string(want));
want = want(~cellfun(@(s) isempty(strtrim(s)), want));
if isempty(want)
    u = unique(cond, 'stable');
    mask = cond == u(1);
    return
end
mask = false(numel(cond), 1);
for i = 1:numel(want)
    pat = strtrim(want{i});
    if any(pat == '*' | pat == '?')
        rx = ['^', regexptranslate('wildcard', pat), '$'];
        mask = mask | ~cellfun('isempty', regexp(cellstr(cond), rx, 'once'));
    else
        mask = mask | (cond == string(pat));
    end
end
end


function lbl = local_cond_label(cond, mask)
u = unique(cond(mask), 'stable');
if isempty(u), lbl = ""; elseif isscalar(u), lbl = u(1); else, lbl = strjoin(u, '+'); end
end


function mag = local_size_metric(scores, sizeby)
switch lower(char(sizeby))
    case 'pos',    mag = max(scores, 0);
    case 'signed', mag = max(scores, 0);   % signed still sizes by positive part
    otherwise,     mag = abs(scores);
end
end


function pos = local_wordcloud_layout(ax, terms, weight, gmax, fr)
% Greedy collision-avoidance packing (largest word first; spiral outward to the
% first slot whose bounding box overlaps nothing already placed), with each
% word measured at its MAXIMUM font size. Produces a tight, non-overlapping
% layout like a static wordcloud; positions are then held fixed for the movie.
k = numel(terms);
weight = weight(:);
maxfs = fr(1) + (fr(2) - fr(1)) * (max(weight, 0) / gmax);
maxfs = max(maxfs, fr(1));

% measure each word's bounding box (data units) at its max font
ht = gobjects(k, 1);
for i = 1:k
    ht(i) = text(ax, 0.5, 0.5, char(terms{i}), 'FontSize', maxfs(i), ...
        'FontWeight', 'bold', 'Interpreter', 'none', 'Visible', 'off');
end
drawnow;
bw = zeros(k, 1); bh = zeros(k, 1);
for i = 1:k
    e = get(ht(i), 'Extent'); bw(i) = e(3); bh(i) = e(4);
end
delete(ht);

pad = 0.010;
[~, ord] = sort(weight, 'descend');       % place biggest first (centre-out)
pos = repmat([0.5 0.5], k, 1);
placed = zeros(0, 4);                       % [cx cy w h]
for oi = 1:k
    i = ord(oi);
    c = local_place_box(placed, bw(i) + pad, bh(i) + pad);
    pos(i, :) = c;
    placed(end + 1, :) = [c, bw(i) + pad, bh(i) + pad]; %#ok<AGROW>
end
end


function c = local_place_box(placed, bw, bh)
% First non-overlapping centre on an outward spiral (inside the canvas if
% possible; relax the canvas bound only if nothing fits inside).
if isempty(placed), c = [0.5 0.5]; return; end
for r = 0:0.006:0.8
    npt = max(12, round(2 * pi * r / 0.008));
    for a = linspace(0, 2 * pi, npt)
        c = [0.5 + r * cos(a), 0.5 + r * sin(a)];
        if c(1) - bw / 2 < 0.02 || c(1) + bw / 2 > 0.98 || ...
                c(2) - bh / 2 < 0.02 || c(2) + bh / 2 > 0.98
            continue
        end
        if ~local_box_overlaps(c, bw, bh, placed), return; end
    end
end
for r = 0:0.01:1.6
    npt = max(12, round(2 * pi * r / 0.01));
    for a = linspace(0, 2 * pi, npt)
        c = [0.5 + r * cos(a), 0.5 + r * sin(a)];
        if ~local_box_overlaps(c, bw, bh, placed), return; end
    end
end
c = [0.5 0.5];
end


function tf = local_box_overlaps(c, bw, bh, placed)
tf = any(abs(placed(:, 1) - c(1)) < (placed(:, 3) + bw) / 2 & ...
         abs(placed(:, 2) - c(2)) < (placed(:, 4) + bh) / 2);
end


function c = local_sign_color(s, clim)
t = min(abs(s) / clim, 1);
if s >= 0
    c = [0.55 + 0.45 * t, 0.25 * (1 - t), 0.20 * (1 - t)];   % red family
else
    c = [0.20 * (1 - t), 0.30 * (1 - t), 0.55 + 0.45 * t];   % blue family
end
end


% ---- video helpers ------------------------------------------------------
function vw = local_open_video(file, fps)
vw = struct('type', 'none', 'obj', [], 'file', char(file), 'fps', fps, 'frames', {{}});
f = char(file);
if isempty(f), return; end
[~, ~, e] = fileparts(f); e = lower(e);
if strcmp(e, '.gif')
    vw.type = 'gif';                              % accumulate; write once at close
elseif strcmp(e, '.avi')
    vw.obj = VideoWriter(f, 'Motion JPEG AVI'); vw.obj.FrameRate = fps; open(vw.obj); vw.type = 'avi';
else
    try
        vw.obj = VideoWriter(f, 'MPEG-4');
    catch
        f = regexprep(f, '\.mp4$', '.avi'); vw.file = f;
        vw.obj = VideoWriter(f, 'Motion JPEG AVI');
    end
    vw.obj.FrameRate = fps; open(vw.obj); vw.type = 'mp4';
end
end

function vw = local_write_frame(vw, fig)
if strcmp(vw.type, 'none'), return; end
frame = getframe(fig);
if strcmp(vw.type, 'gif')
    vw.frames{end + 1} = frame2im(frame);        % keep truecolor; quantize at close
else
    writeVideo(vw.obj, frame);
end
end

function vw = local_close_video(vw)
if strcmp(vw.type, 'gif')
    F = vw.frames;
    if isempty(F), return; end
    % ONE global colormap from all frames, so colors (red/blue/grey) survive
    % even when the first frame is near-blank (no significant terms yet).
    [~, gmap] = rgb2ind(cat(1, F{:}), 256);
    dt = 1 / max(vw.fps, 1);
    for i = 1:numel(F)
        A = rgb2ind(F{i}, gmap);
        if i == 1
            imwrite(A, gmap, vw.file, 'gif', 'LoopCount', Inf, 'DelayTime', dt);
        else
            imwrite(A, gmap, vw.file, 'gif', 'WriteMode', 'append', 'DelayTime', dt);
        end
    end
elseif ~strcmp(vw.type, 'none') && ~isempty(vw.obj)
    close(vw.obj);
end
end

function pth = local_video_path(vw, file)
if isempty(char(file)), pth = ''; else, pth = vw.file; end
end

function s = local_file_note(f)
if isempty(f), s = ' (not saved)'; else, s = sprintf(' -> %s', f); end
end
