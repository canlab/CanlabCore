function fig = hrf_plot_causality(R, varargin)
%HRF_PLOT_CAUSALITY Visualize a directed (Granger) net-flow result.
%
% Three complementary views of an hrf_causality / hrf_causality_analyze
% result, in one figure:
%   A  net-flow heatmap        - net(i,j) = i->j minus j->i (antisymmetric),
%                                diverging colormap, significant cells boxed.
%   B  directed graph          - circular layout; an arrow i->j for each
%                                significant POSITIVE net edge, width/colour
%                                by magnitude (the "who leads whom" picture).
%   C  net-drive ranking       - per-node net outflow (row sum of net):
%                                positive = net source/driver, negative = sink.
%
% The net matrix is antisymmetric, so the heatmap's two triangles mirror
% (sign-flipped) and the graph only needs the positive direction.
%
% :Usage:
% ::
%     hrf_plot_causality(R)                       % first evoked mode
%     hrf_plot_causality(R, 'Mode','remove', 'PThresh',0.05)
%     hrf_plot_causality(R, 'PField','p', 'TopEdges',30, 'Save','flow.png')
%     hrf_plot_causality(net_matrix, 'Nodes', names)   % bare matrix
%
% :Inputs:
%   **R:** an hrf_causality(_analyze) struct (uses R.(Mode).net_group / .p_fdr
%          / .p and R.nodes), a struct with fields net_group + p/p_fdr +
%          nodes, or a bare [N x N] net matrix.
%
% :Optional Inputs:
%   **'Mode':**    evoked mode field to plot ('remove'/'keep'); default = first
%                  in R.modes.
%   **'PThresh':** significance threshold for edges/boxes. Default 0.05.
%   **'PField':**  'p_fdr' (default) or 'p'. Ignored if absent.
%   **'Nodes':**   node names (for a bare-matrix input or to override).
%   **'TopEdges':**cap the graph to the strongest N significant edges. If
%                  NONE are significant, the strongest N by |net| are drawn
%                  dashed (with a note). Default 40.
%   **'Title':**   figure title. Default 'HRF causality (net directed flow)'.
%   **'Save':**    path to write a PNG. Default ''.
%
% :Output:  **fig** - the figure handle.
%
% See also: hrf_causality, hrf_causality_analyze, hrf_granger_causality.

p = inputParser;
p.addRequired('R');
p.addParameter('Mode', '', @(x) ischar(x) || isstring(x));
p.addParameter('PThresh', 0.05, @(x) isscalar(x) && x > 0);
p.addParameter('PField', 'p_fdr', @(x) ischar(x) || isstring(x));
p.addParameter('Nodes', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('TopEdges', 40, @(x) isscalar(x) && x >= 1);
p.addParameter('Title', 'HRF causality (net directed flow)', @(x) ischar(x) || isstring(x));
p.addParameter('Save', '', @(x) ischar(x) || isstring(x));
p.parse(R, varargin{:});
opts = p.Results;

[net, pv, nodes, modelabel] = local_unpack(R, opts);
N = size(net, 1);
net(1:N+1:end) = 0;
if isempty(pv), pv = nan(N); end
sig = pv < opts.PThresh;
sig(1:N+1:end) = false;

fig = figure('Color', 'w', 'Position', [60 60 1500 560], 'Name', char(opts.Title));
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ttl = char(opts.Title); if ~isempty(modelabel), ttl = sprintf('%s  [%s]', ttl, modelabel); end
title(tl, ttl, 'FontWeight', 'bold', 'Interpreter', 'none');

% ---- A: heatmap --------------------------------------------------------
ax = nexttile(tl);
imagesc(ax, net); axis(ax, 'square');
colormap(ax, local_diverging(256));
lim = max(abs(net(:))); if lim == 0 || ~isfinite(lim), lim = 1; end
set(ax, 'CLim', [-lim lim]);
cb = colorbar(ax); cb.Label.String = 'net flow  (row \rightarrow col)';
set(ax, 'XTick', 1:N, 'XTickLabel', nodes, 'YTick', 1:N, 'YTickLabel', nodes, ...
    'TickLabelInterpreter', 'none', 'FontSize', 7, 'XTickLabelRotation', 90);
hold(ax, 'on');
[ri, ci] = find(triu(sig, 1) | tril(sig, -1));
plot(ax, ci, ri, 's', 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 8, 'LineWidth', 0.75);
title(ax, sprintf('net-flow matrix (boxed: %s < %.2g)', char(opts.PField), opts.PThresh), ...
    'FontSize', 9, 'Interpreter', 'none');

% ---- B: directed graph -------------------------------------------------
ax = nexttile(tl);
local_flow_graph(ax, net, sig, nodes, opts.TopEdges);

% ---- C: net-drive ranking ----------------------------------------------
ax = nexttile(tl);
drive = sum(net, 2);                          % net outflow per node
[ds, ord] = sort(drive, 'ascend');
b = barh(ax, ds, 'FaceColor', 'flat');
cmap = local_diverging(256);
b.CData = local_map_colors(ds, cmap);
set(ax, 'YTick', 1:N, 'YTickLabel', nodes(ord), 'TickLabelInterpreter', 'none', 'FontSize', 7);
ylim(ax, [0.5 N+0.5]); grid(ax, 'on'); box(ax, 'off');
xline(ax, 0, 'Color', [0.4 0.4 0.4]);
xlabel(ax, 'net outflow  (\Sigma_j net_{ij})', 'Interpreter', 'tex');
title(ax, 'net driver (+)  \leftrightarrow  receiver (−)', 'FontSize', 9, 'Interpreter', 'tex');

if ~isempty(char(opts.Save))
    try, exportgraphics(fig, char(opts.Save), 'Resolution', 130); catch, saveas(fig, char(opts.Save)); end
end
end


% =========================================================================
function local_flow_graph(ax, net, sig, nodes, topEdges)
N = size(net, 1);
mask = sig & (net > 0);
dashed = false;
if ~any(mask(:))
    % nothing significant: fall back to strongest positive edges, dashed
    dashed = true;
    pos = net; pos(pos <= 0) = 0;
    thr = local_kth_largest(pos(:), topEdges);
    mask = pos >= thr & pos > 0;
end
[si, ti] = find(mask);
w = net(sub2ind([N N], si, ti));
% cap to strongest topEdges
if numel(w) > topEdges
    [~, k] = sort(w, 'descend'); k = k(1:topEdges);
    si = si(k); ti = ti(k); w = w(k);
end
if isempty(si)
    axis(ax, 'off'); text(ax, 0.5, 0.5, 'no edges to draw', 'HorizontalAlignment', 'center');
    return
end
G = digraph(si, ti, w, N);
ang = linspace(0, 2*pi, N+1); ang = ang(1:N);
xy = [cos(ang)', sin(ang)'];
lw = 0.5 + 4.5 * (abs(w) / max(abs(w)));
h = plot(ax, G, 'XData', xy(:,1), 'YData', xy(:,2), 'NodeLabel', nodes, ...
    'LineWidth', lw, 'ArrowSize', 9, 'Interpreter', 'none', 'NodeColor', [0.2 0.2 0.2], 'NodeFontSize', 7);
% colour edges by weight
ew = G.Edges.Weight;
cmap = local_diverging(256);
h.EdgeColor = local_map_colors(ew, cmap);
if dashed, h.LineStyle = '--'; end
axis(ax, 'equal', 'off');
note = ''; if dashed, note = sprintf('  (none sig; top %d by |net|)', topEdges); end
title(ax, sprintf('directed flow  i\\rightarrowj%s', note), 'FontSize', 9, 'Interpreter', 'tex');
end


function [net, pv, nodes, modelabel] = local_unpack(R, opts)
modelabel = '';
pv = [];
if isnumeric(R)
    net = R; nodes = local_default_nodes(opts.Nodes, size(R, 1));
    return
end
if ~isstruct(R), error('hrf_plot_causality:Input', 'R must be a struct or a matrix.'); end
if isfield(R, 'modes')
    mode = char(opts.Mode); if isempty(mode), mode = R.modes{1}; end
    if ~isfield(R, mode), error('hrf_plot_causality:Mode', 'R has no mode ''%s''.', mode); end
    S = R.(mode); modelabel = mode;
elseif isfield(R, 'net_group')
    S = R;
else
    error('hrf_plot_causality:Struct', 'Struct must have .modes or .net_group.');
end
net = S.net_group;
pf = char(opts.PField);
if isfield(S, pf), pv = S.(pf); elseif isfield(S, 'p'), pv = S.p; end
if ~isempty(opts.Nodes)
    nodes = cellstr(string(opts.Nodes));
elseif isfield(R, 'nodes') && ~isempty(R.nodes)
    nodes = cellstr(string(R.nodes));
elseif isfield(S, 'nodes') && ~isempty(S.nodes)
    nodes = cellstr(string(S.nodes));
else
    nodes = local_default_nodes({}, size(net, 1));
end
end


function nodes = local_default_nodes(want, N)
if ~isempty(want), nodes = cellstr(string(want)); else, nodes = arrayfun(@(i) sprintf('n%d', i), 1:N, 'uni', 0); end
end

function C = local_diverging(n)
% blue -> white -> red
half = floor(n/2);
up = linspace(0, 1, half)';
top = [ones(n-half,1), linspace(1,0,n-half)', linspace(1,0,n-half)'];
bot = [up, up, ones(half,1)];
C = [bot; top];
end

function C = local_map_colors(v, cmap)
v = v(:); lim = max(abs(v)); if lim == 0 || ~isfinite(lim), lim = 1; end
idx = round((v + lim) / (2*lim) * (size(cmap,1)-1)) + 1;
idx = min(max(idx, 1), size(cmap,1));
C = cmap(idx, :);
end

function thr = local_kth_largest(x, k)
x = sort(x(:), 'descend');
k = min(k, numel(x));
thr = x(max(k,1));
end
