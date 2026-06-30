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
%   **'TopN':**       number of terms shown (ranked by peak |score| over lags).
%                     Default 35.
%   **'SizeBy':**     'abs' (default) | 'pos' | 'signed' -- which magnitude
%                     drives font size (abs = both directions grow).
%   **'FrameRate':**  fps of the movie. Default 4.
%   **'OutputFile':** '.mp4'/'.avi'/'.gif' to write; '' = just show. Default ''.
%   **'FontRange':**  [min max] point sizes. Default [9 46].
%   **'Title':**      title prefix. Default 'Neurosynth terms over the HRF'.
%   **'Verbose'/'doverbose':** chatter. Default true.
%
% :Outputs:
%   **out:** struct with .scores [nLag x nTerm_kept], .terms, .lags, .pos
%            (fixed [k x 2] layout), .file (written path or '').
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
p.addParameter('TopN', 35, @(x) isscalar(x) && x >= 1);
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

[scores, terms, lags, condlabel] = local_get_term_matrix(source, opts);   % scores [nLag x nTerm]
if isempty(scores), error('hrf_animate_wordcloud:NoData', 'No map_%s_* scores found for the requested condition.', char(opts.Set)); end

% rank + keep top terms by peak magnitude over lags
mag = local_size_metric(scores, opts.SizeBy);
imp = max(mag, [], 1, 'omitnan');
[~, ord] = sort(imp, 'descend');
keep = ord(1:min(opts.TopN, numel(terms)));
terms = terms(keep); scores = scores(:, keep); mag = mag(:, keep);
pos = local_spiral_positions(numel(terms));        % fixed layout, importance-ordered

gmax = max(mag(:)); if gmax == 0 || ~isfinite(gmax), gmax = 1; end
clim = max(abs(scores(:))); if clim == 0 || ~isfinite(clim), clim = 1; end
fr = opts.FontRange;

fig = figure('Color', 'w', 'Position', [80 80 900 700], 'Name', char(opts.Title));
ax = axes(fig, 'Position', [0.02 0.05 0.96 0.88]); axis(ax, [0 1 0 1]); axis(ax, 'off'); hold(ax, 'on');

vw = local_open_video(opts.OutputFile, opts.FrameRate);
nLag = numel(lags);
for li = 1:nLag
    cla(ax); axis(ax, [0 1 0 1]); axis(ax, 'off');
    for ti = 1:numel(terms)
        s = scores(li, ti); m = mag(li, ti);
        if ~isfinite(m) || m <= 0, continue; end
        fs = fr(1) + (fr(2) - fr(1)) * (m / gmax);
        text(ax, pos(ti, 1), pos(ti, 2), char(terms{ti}), 'FontSize', fs, ...
            'Color', local_sign_color(s, clim), 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'Interpreter', 'none', 'Clipping', 'on');
    end
    title(ax, sprintf('%s   —   %s   t = %.1f s', char(opts.Title), condlabel, lags(li)), ...
        'Interpreter', 'none', 'FontSize', 12);
    drawnow;
    vw = local_write_frame(vw, fig);
end
vw = local_close_video(vw);

out = struct('scores', scores, 'terms', {terms}, 'lags', lags(:)', 'pos', pos, ...
    'file', local_video_path(vw, opts.OutputFile));
if verbose
    fprintf('hrf_animate_wordcloud: %d terms x %d lags (condition %s)%s\n', ...
        numel(terms), nLag, condlabel, local_file_note(out.file));
end
end


% =========================================================================
function [scores, terms, lags, condlabel] = local_get_term_matrix(source, opts)
prefix = ['map_', char(opts.Set), '_'];
if isstruct(source) && isfield(source, 'scores')
    scores = source.scores; terms = cellstr(string(source.terms));
    lags = source.lags(:)'; condlabel = 'curve';
    return
end
% cell of dirs / input tables, or a single output DIRECTORY -> pool into one
% input table (group-mean across every subject of every dir), then proceed.
is_dir = (ischar(source) || isstring(source)) && isfolder(char(string(source)));
if iscell(source) || is_dir
    IT = local_pool_input_table(source);
    [scores, terms, lags, condlabel] = local_matrix_from_input_table(IT, prefix, opts);
    return
end
if (ischar(source) || isstring(source)) && endsWith(string(source), '.csv')
    [scores, terms, lags, condlabel] = local_matrix_from_csv(char(source), prefix, opts);
    return
end
if istable(source)
    [scores, terms, lags, condlabel] = local_matrix_from_input_table(source, prefix, opts);
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


function [M, terms, lags, condlabel] = local_matrix_from_input_table(IT, prefix, opts)
% Pool group-mean across subjects: average each (term, lag) over rows' CSVs.
file_col = '';
if strcmpi(char(opts.Object), 't') && any(strcmp('t_scores_file', IT.Properties.VariableNames))
    file_col = 't_scores_file';
elseif any(strcmp('beta_scores_file', IT.Properties.VariableNames))
    file_col = 'beta_scores_file';
end
if isempty(file_col), error('hrf_animate_wordcloud:NoScoreFiles', 'input_table lacks *_scores_file columns.'); end
model = char(opts.Model);   % strcmpi below is case-insensitive
acc = []; terms = {}; lags = []; condlabel = '';
nfile = 0;
for i = 1:height(IT)
    if any(strcmp('model', IT.Properties.VariableNames))
        if ~strcmpi(char(string(IT.model(i))), model), continue; end
    end
    f = char(string(IT.(file_col)(i)));
    if isempty(f) || exist(f, 'file') ~= 2, continue; end
    T = readtable(f, 'TextType', 'string');
    [Mi, ti, li, cl] = local_matrix_from_table(T, prefix, opts);
    if isempty(Mi), continue; end
    if isempty(acc)
        acc = Mi; terms = ti; lags = li; condlabel = cl; nfile = 1;
    elseif isequal(ti, terms) && isequal(li, lags)
        acc = acc + Mi; nfile = nfile + 1;
    else
        warning('hrf_animate_wordcloud:GridMismatch', ...
            'Skipping a score file whose term/lag grid differs from the first.');
    end
end
if isempty(acc), M = []; return; end
M = acc / max(nfile, 1);
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


function pos = local_spiral_positions(k)
% Archimedean spiral in [0,1]^2, densest (most important) near centre.
pos = zeros(k, 2);
golden = pi * (3 - sqrt(5));
for i = 1:k
    r = 0.46 * sqrt(i / max(k, 1));
    th = (i - 1) * golden;
    pos(i, :) = [0.5 + r * cos(th), 0.5 + r * sin(th)];
end
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
vw = struct('type', 'none', 'obj', [], 'file', char(file), 'gifinit', false);
f = char(file);
if isempty(f), return; end
[~, ~, e] = fileparts(f); e = lower(e);
if strcmp(e, '.gif')
    vw.type = 'gif';
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
    [A, map] = rgb2ind(frame2im(frame), 256);
    if ~vw.gifinit
        imwrite(A, map, vw.file, 'gif', 'LoopCount', Inf, 'DelayTime', 0.25); vw.gifinit = true;
    else
        imwrite(A, map, vw.file, 'gif', 'WriteMode', 'append', 'DelayTime', 0.25);
    end
else
    writeVideo(vw.obj, frame);
end
end

function vw = local_close_video(vw)
if ~strcmp(vw.type, 'none') && ~strcmp(vw.type, 'gif') && ~isempty(vw.obj)
    close(vw.obj);
end
end

function pth = local_video_path(vw, file)
if isempty(char(file)), pth = ''; else, pth = vw.file; end
end

function s = local_file_note(f)
if isempty(f), s = ' (not saved)'; else, s = sprintf(' -> %s', f); end
end
