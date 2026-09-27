function h = plot(obj, varargin)
% plot  Visualize an rsm.
%
% Modes
%   plot(R)                    % raw heatmap (default) — actual data values
%   plot(R, 'mode','rank')     % rank-transformed heatmap (Kriegeskorte style)
%   plot(R, 'mode','mds')      % 2D MDS scatter of conditions
%   plot(R, 'mode','dendrogram') % hierarchical clustering on rows
%   plot(R, 'mode','grid')     % subplot grid, one per replicate
%
%   plot(R, 'subject', 3)        % plot the 3rd slice along dim 3
%   plot(R, 'matched_pairs', P)  % overlay white rectangles around (i,j) pairs
%                                % (P must have indices in 1..k; out-of-range warns + drops)
%   plot(R, 'block_borders_by', 'bodysite')   % borders around metadata blocks
%   plot(R, 'border_color', 'white')   % 'white' (default), 'black', 'auto', RGB, or Nx3
%   plot(R, 'climits', [-1 1])   % override data-driven color limits
%   plot(R, 'colormap', cmap)    % override default colormap
%   plot(R, 'title', 'My RSM')   % override title
%
% Default behavior (mode='heatmap')
%   - RSMs (similarity matrices, signed values): diverging blue-white-red
%     colormap (colormap_tor) with symmetric color limits [-max|x|, +max|x|].
%   - RDMs (dissimilarity matrices, non-negative values): sequential parula
%     colormap with limits [0, max(x)].
%
% Returns a struct of handles (cb = colorbar, ax = axes, ...).

if builtin('numel', obj) > 1
    % Array of rsm objects (e.g., parcelwise output). Auto-grid each entry
    % into subplots. Caps at 16 entries to keep things readable; for larger
    % arrays the user should slice or use get_by_label / select.
    h = plot_array(obj, varargin{:});
    return
end

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('mode',              'heatmap', @(x) ischar(x) || isstring(x));
p.addParameter('subject',           [],        @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('matched_pairs',     [],        @(x) isempty(x) || (isnumeric(x) && size(x,2)==2));
p.addParameter('block_borders_by',  '',        @(x) ischar(x) || isstring(x));
p.addParameter('border_color',      'white',   @(x) ischar(x) || isstring(x) || (isnumeric(x) && (numel(x)==3 || size(x,2)==3)));
p.addParameter('border_width',      2.5,       @isnumeric);
p.addParameter('climits',           [],        @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
p.addParameter('colormap',          [],        @(x) isempty(x) || (isnumeric(x) && size(x,2)==3) || ischar(x) || isstring(x));
p.addParameter('title',             '',        @(x) ischar(x) || isstring(x));
p.addParameter('label_fontsize',    10,        @isnumeric);
p.addParameter('rotate_xticks',     true,      @(x) islogical(x) || isnumeric(x));
p.addParameter('mds_engine',        'builtin', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

mode = lower(char(p.Results.mode));

switch mode
    case 'grid'
        h = plot_grid(obj, p.Results);
    case 'mds'
        h = plot_mds(obj, p.Results);
    case 'dendrogram'
        h = plot_dendrogram(obj, p.Results);
    case 'heatmap'
        h = plot_heatmap(obj, p.Results, 'raw');
    case 'raw'   % legacy alias
        h = plot_heatmap(obj, p.Results, 'raw');
    case 'rank'
        h = plot_heatmap(obj, p.Results, 'rank');
    otherwise
        error('rsm:plot:badMode', ['Unknown mode: %s. ', ...
            'Valid: heatmap (default), rank, mds, dendrogram, grid.'], mode);
end

end


% =========================================================================
function h = plot_array(obj, varargin)

n = builtin('numel', obj);
max_show = 16;

if n > max_show
    warning('rsm:plot:tooManyArrayEntries', ...
        ['plot() on an %d-element rsm array shows the first %d only. ', ...
         'Use obj.get_by_label(name), obj(idx), or arrayfun(@(R) plot(R), obj) for full control.'], ...
        n, max_show);
    n_show = max_show;
else
    n_show = n;
end

ncols = ceil(sqrt(n_show)); nrows = ceil(n_show / ncols);
h.fig = gcf;
h.subaxes = gobjects(n_show, 1);

for i = 1:n_show
    h.subaxes(i) = subplot(nrows, ncols, i);
    % Pull each rsm's source as the subplot title
    plot(obj(i), 'title', short_source(obj(i)), varargin{:});
end

end


function s = short_source(R)
% Pick a concise title from an rsm's source field
s = char(R.source);
if startsWith(s, 'parcel:'), s = s(8:end); end   % strip prefix
if numel(s) > 28, s = [s(1:25) '...']; end
% Use CanlabCore's label prettifier when available
if exist('format_strings_for_legend', 'file') == 2
    try, s = format_strings_for_legend({s}); s = s{1}; catch, end
end
end


function pretty = prettify_labels(labels)
% Run labels through format_strings_for_legend (CanlabCore utility) so
% underscores render as readable subscripts in MATLAB plot ticks. Falls
% back to raw labels if the function isn't on the path.

if exist('format_strings_for_legend', 'file') ~= 2
    pretty = labels;
    return
end
try
    pretty = format_strings_for_legend(labels);
catch
    pretty = labels;
end
end


% =========================================================================
function h = plot_heatmap(obj, opts, render_mode)
% render_mode  'raw'  -- show actual data values (default)
%              'rank' -- rank-transformed for Kriegeskorte-style display

% Slice or aggregate down to a 2D matrix
M = collapse_to_2d(obj, opts.subject);

caps = probe_rsatoolbox();

is_rank = strcmp(render_mode, 'rank');

if is_rank
    if caps.rdm_rankTransform && caps.rdm_squareRDM
        try
            if obj.is_dissimilarity
                M_disp = rankTransform_safe(M, caps);
            else
                rdm = 1 - M;
                M_disp = rankTransform_safe(rdm, caps);
            end
        catch
            M_disp = local_rank_transform(M);
        end
    else
        M_disp = local_rank_transform(M);
    end
    if isempty(opts.climits), clim = [0 1]; else, clim = opts.climits; end
else
    % Raw mode -- show actual values; auto-pick climits and colormap
    M_disp = M;
    [clim, default_cmap_kind] = autoscale_for_raw(M_disp, obj.is_dissimilarity, opts.climits);
end

% Render
h.fig = gcf;
h.ax  = gca;
h.im  = imagesc(M_disp, clim);
axis image;
h.cb  = colorbar;

% Colormap selection
if ~isempty(opts.colormap)
    if isnumeric(opts.colormap)
        cmap = opts.colormap;
    else
        cmap = feval(char(opts.colormap), 256);
    end
elseif is_rank
    % Rank display uses sequential colormap
    if caps.util_RDMcolormap
        try, cmap = RDMcolormap; catch, cmap = parula(256); end
    else
        cmap = parula(256);
    end
else
    % Raw display: diverging for signed (RSM), sequential for non-negative (RDM)
    cmap = pick_default_cmap(default_cmap_kind);
end
colormap(h.ax, cmap);

% Labels
k = size(M_disp, 1);
if ~isempty(obj.labels)
    pretty_labels = prettify_labels(obj.labels);
    set(h.ax, 'XTick', 1:k, 'XTickLabel', pretty_labels, ...
              'YTick', 1:k, 'YTickLabel', pretty_labels, ...
              'FontSize', opts.label_fontsize);
    if opts.rotate_xticks, set(h.ax, 'XTickLabelRotation', 45); end
end

% Title
ttl = char(opts.title);
if isempty(ttl)
    kind = 'RSM'; if obj.is_dissimilarity, kind = 'RDM'; end
    ttl = sprintf('%s (%s%s)', kind, obj.metric, ...
        cond_subtitle(obj, opts.subject));
end
title(h.ax, ttl, 'Interpreter', 'none');

% Optional overlays
if ~isempty(opts.matched_pairs)
    P = opts.matched_pairs;
    % Validate: drop any out-of-range pairs with a warning
    k = size(M_disp, 1);
    bad = any(P < 1, 2) | any(P > k, 2);
    if any(bad)
        warning('rsm:plot:matchedPairsOutOfRange', ...
            ['matched_pairs contained %d pair(s) with indices outside 1..%d ', ...
             '(rsm is %dx%d); those will be dropped.'], sum(bad), k, k, k);
        P = P(~bad, :);
    end
    if ~isempty(P)
        hold(h.ax, 'on');
        for ii = 1:size(P, 1)
            r = P(ii, 1); c = P(ii, 2);
            rectangle('Position', [c-0.5, r-0.5, 1, 1], 'EdgeColor', 'w', 'LineWidth', 2);
            rectangle('Position', [r-0.5, c-0.5, 1, 1], 'EdgeColor', 'w', 'LineWidth', 2);
        end
        hold(h.ax, 'off');
    end
end

if ~isempty(opts.block_borders_by)
    overlay_block_borders(h.ax, obj, char(opts.block_borders_by), ...
        opts.border_color, opts.border_width);
end

end


% =========================================================================
function h = plot_grid(obj, opts) %#ok<INUSD>

if size(obj.dat, 3) <= 1
    error('rsm:plot:gridNonReplicate', 'mode=''grid'' requires a 3D rsm (size(dat,3) > 1).');
end

N = size(obj.dat, 3);
ncols = ceil(sqrt(N)); nrows = ceil(N / ncols);
h.fig = gcf;
h.subaxes = gobjects(N, 1);

for n = 1:N
    h.subaxes(n) = subplot(nrows, ncols, n);
    slice = obj.dat(:, :, n);
    slice_obj = rsm(slice, 'metric', obj.metric, 'labels', obj.labels, ...
                    'is_dissimilarity', obj.is_dissimilarity, ...
                    'groupings', obj.groupings);
    % Recurse on the single slice
    plot(slice_obj, 'mode', 'heatmap', 'title', replicate_title(obj, n));
end

end


% =========================================================================
function h = plot_mds(obj, opts)

M = collapse_to_2d(obj, []);
if ~obj.is_dissimilarity, M = 1 - M; end

% Engine: 'builtin' (default) gives a clean labeled cmdscale scatter that
% reads well at any condition count. 'rsatoolbox' routes through
% rsa.fig.MDSConditions when available (its multi-panel display is busier,
% especially for few conditions).
engine = 'builtin';
if isfield(opts, 'mds_engine') && ~isempty(opts.mds_engine)
    engine = lower(char(opts.mds_engine));
end

if strcmp(engine, 'rsatoolbox')
    caps = probe_rsatoolbox();
    if caps.fig_MDSConditions && ~isempty(which('rsa.fig.MDSConditions'))
        try
            rsa.fig.MDSConditions(struct('RDM', M, 'name', 'rsm'), struct('saveFiguresFig', false));
            h.fig = gcf;
            return
        catch
            % fall through to builtin
        end
    end
end

% Built-in classical MDS via cmdscale (default)
M = (M + M') / 2;
M(1:size(M,1)+1:end) = 0;
[Y, e] = cmdscale(M);
if size(Y, 2) < 2, Y = [Y, zeros(size(Y,1), 2 - size(Y,2))]; end
Y = Y(:, 1:2);

h.fig = gcf;
h.ax  = gca;
h.scatter = scatter(Y(:,1), Y(:,2), 90, 'filled');
hold(h.ax, 'on');
if ~isempty(obj.labels)
    pretty = prettify_labels(obj.labels);
    for i = 1:size(Y, 1)
        text(Y(i,1), Y(i,2), ['  ' pretty{i}], 'FontSize', 10, 'Interpreter', 'none');
    end
end
hold(h.ax, 'off');
axis equal; grid on;
xlabel(h.ax, 'MDS dim 1'); ylabel(h.ax, 'MDS dim 2');
if numel(e) >= 2 && sum(abs(e)) > 0
    var2d = 100 * sum(e(1:2)) / sum(abs(e));
    title(h.ax, sprintf('MDS (%s, %.0f%% in 2D)', obj.metric, var2d), 'Interpreter', 'none');
else
    title(h.ax, sprintf('MDS (%s)', obj.metric), 'Interpreter', 'none');
end

end


% =========================================================================
function h = plot_dendrogram(obj, opts) %#ok<INUSD>

M = collapse_to_2d(obj, []);
if ~obj.is_dissimilarity, M = 1 - M; end

% squareform expects a symmetric matrix with zeros on the diagonal
M_clean = M;
k = size(M_clean, 1);
M_clean(1:k+1:end) = 0;
% Symmetrize to remove any tiny asymmetry
M_clean = (M_clean + M_clean') / 2;
% Floor at 0 in case rounding gave tiny negatives
M_clean(M_clean < 0) = 0;

d = squareform(M_clean, 'tovector');
Z = linkage(d, 'average');

h.fig = gcf; h.ax = gca;
if ~isempty(obj.labels)
    [h.handle, ~, h.order] = dendrogram(Z, 0, 'Labels', obj.labels);
else
    [h.handle, ~, h.order] = dendrogram(Z, 0);
end
title(h.ax, sprintf('Dendrogram (%s, average linkage)', obj.metric), 'Interpreter', 'none');

end


% =========================================================================
function M = collapse_to_2d(obj, subject_idx)
% Pick a 2D slice from a possibly-3D rsm.

if size(obj.dat, 3) == 1
    M = obj.dat;
elseif ~isempty(subject_idx)
    M = obj.dat(:, :, subject_idx);
else
    M = mean(obj.dat, 3, 'omitnan');
end

end


function s = cond_subtitle(obj, subject_idx)
if size(obj.dat, 3) == 1
    s = '';
elseif ~isempty(subject_idx)
    s = sprintf(', slice %d', subject_idx);
else
    s = ', mean across replicates';
end
end


function s = replicate_title(obj, n)
if ~isempty(obj.replicate_table) && n <= height(obj.replicate_table)
    parts = {};
    for v = obj.replicate_table.Properties.VariableNames
        val = obj.replicate_table.(v{1})(n);
        parts{end+1} = sprintf('%s=%s', v{1}, char(string(val))); %#ok<AGROW>
    end
    s = strjoin(parts, ', ');
else
    s = sprintf('replicate %d', n);
end
end


function [clim, cmap_kind] = autoscale_for_raw(M, is_dissim, user_clim)
% Data-driven color limits + colormap kind for raw heatmap display.
%
% cmap_kind = 'diverging'  -> blue-white-red (signed values, RSMs)
%             'sequential' -> parula (non-negative values, RDMs / distances)

vals = M(:);
vals = vals(isfinite(vals));

if isempty(vals)
    clim = [0 1]; cmap_kind = 'sequential'; return
end

if ~isempty(user_clim)
    clim = user_clim;
    % Choose colormap from the requested limits
    if clim(1) < 0 && clim(2) > 0, cmap_kind = 'diverging';
    else, cmap_kind = 'sequential';
    end
    return
end

vmin = min(vals); vmax = max(vals);

if is_dissim || vmin >= -1e-9
    % Non-negative distances or essentially-non-negative RSM
    clim = [max(0, vmin), vmax];
    if clim(1) == clim(2), clim(2) = clim(1) + 1; end
    cmap_kind = 'sequential';
else
    % Signed values -> symmetric limits
    amax = max(abs([vmin vmax]));
    clim = [-amax, amax];
    cmap_kind = 'diverging';
end

end


function cmap = pick_default_cmap(kind)
% Default colormap selector. For 'diverging' tries colormap_tor (CanlabCore
% blue-white-red) and falls back to a built-in blue-white-red.

switch kind
    case 'diverging'
        if exist('colormap_tor', 'file') == 2
            try
                cmap = colormap_tor([0 0 1], [1 0 0], [1 1 1]);
                if size(cmap, 1) < 64, cmap = interp_cmap(cmap, 256); end
                return
            catch
                % fall through
            end
        end
        % Stock fallback: 3-stop bwr
        n = 256;
        half = floor(n/2);
        b2w = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
        w2r = [ones(n-half,1), linspace(1,0,n-half)', linspace(1,0,n-half)'];
        cmap = [b2w; w2r];
    case 'sequential'
        cmap = parula(256);
    otherwise
        cmap = parula(256);
end

end


function cmap_out = interp_cmap(cmap_in, n_out)
% Linearly interpolate a colormap to n_out rows.
n_in = size(cmap_in, 1);
xi = linspace(1, n_in, n_out);
cmap_out = zeros(n_out, 3);
for c = 1:3, cmap_out(:, c) = interp1(1:n_in, cmap_in(:, c), xi)'; end
end


function M_t = rankTransform_safe(M, caps)
if ~isempty(which('rsa.rdm.rankTransform')) && ~isempty(which('rsa.rdm.squareRDM'))
    M_t = rsa.rdm.rankTransform(rsa.rdm.squareRDM(M), 1);
elseif caps.rdm_rankTransform && caps.rdm_squareRDM
    M_t = rankTransform(squareRDM(M), 1);
else
    M_t = local_rank_transform(M);
end
n = size(M_t, 1);
M_t(1:n+1:end) = 1;
end


function overlay_block_borders(ax, obj, col, border_color, border_width)
% Draw rectangles around contiguous blocks of equal values in
% obj.metadata_table.(col).
%
% border_color
%   'white' / 'black' (default 'white') -- single contrasting color
%   'auto'          -- per-block colors from lines() (legacy behavior; can clash with heatmap)
%   [r g b]         -- single RGB triplet for all blocks
%   [N x 3] matrix  -- per-block RGB rows (N == number of blocks)
%
% border_width  scalar line width (default 2.5).

if isempty(obj.metadata_table) || ~ismember(col, obj.metadata_table.Properties.VariableNames)
    warning('rsm:plot:blockBordersMissing', ...
        'metadata_table does not contain column "%s"; skipping block borders.', col);
    return
end

if nargin < 4 || isempty(border_color), border_color = 'white'; end
if nargin < 5 || isempty(border_width), border_width = 2.5; end

v = obj.metadata_table.(col);
[grp, ~] = findgroups(v);
n_blocks = max(grp);

cmap = resolve_border_cmap(border_color, n_blocks);

hold(ax, 'on');
for g = 1:n_blocks
    idx = find(grp == g);
    if isempty(idx), continue; end
    s = min(idx); e = max(idx);
    rectangle('Parent', ax, 'Position', [s-0.5, s-0.5, e-s+1, e-s+1], ...
        'EdgeColor', cmap(g, :), 'LineWidth', border_width);
end
hold(ax, 'off');

end


function cmap = resolve_border_cmap(border_color, n_blocks)
% Returns an [n_blocks x 3] colormap, with the same color repeated when a
% single color is requested.

if ischar(border_color) || isstring(border_color)
    name = lower(char(border_color));
    switch name
        case 'white', single = [1 1 1];
        case 'black', single = [0 0 0];
        case 'auto',  cmap = lines(n_blocks); return
        otherwise
            % Try MATLAB color shorthand ('r','g','b',...)
            try
                single = validatecolor(name);
            catch
                warning('rsm:plot:badBorderColor', ...
                    'Unknown border_color "%s"; falling back to white.', name);
                single = [1 1 1];
            end
    end
    cmap = repmat(single, n_blocks, 1);
elseif isnumeric(border_color)
    if numel(border_color) == 3
        cmap = repmat(border_color(:)', n_blocks, 1);
    elseif size(border_color, 2) == 3 && size(border_color, 1) >= n_blocks
        cmap = border_color(1:n_blocks, :);
    else
        warning('rsm:plot:badBorderColor', ...
            'Numeric border_color must be a 1x3 RGB or Nx3 matrix; using white.');
        cmap = repmat([1 1 1], n_blocks, 1);
    end
end

end
