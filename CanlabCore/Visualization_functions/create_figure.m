function [f1, axh, tag] = create_figure(tagname, varargin)
% Create or reuse a tagged figure, lay out subplots, size it so axes are
% equally sized across calls, and place it without overlapping existing
% figure windows.
%
% :Usage:
% ::
%
%    [f1, axh, tag] = create_figure(['tagname'], [subplotrows], [subplotcols], ...
%                                   [keep_existing_axes], [force_resize])
%
% Looks for an existing figure with the given Tag. If found, reuses (and
% optionally clears) it; otherwise creates a new figure. When creating a
% new figure (or when forced), the figure is sized so each subplot occupies
% roughly the same on-screen area as in any other call -- figures with more
% subplots are larger overall, capped at 80% of the screen in either
% dimension. Newly created figures are placed at the first top-left grid
% position that does not overlap any existing visible figure; if no free
% spot exists, the new figure cascades from the top-left corner.
%
% :Inputs:
%
%   **tagname:**
%        String/char Tag to identify the figure. Default: 'nmdsfig'. If a
%        figure with this Tag exists, it is reused.
%
% :Optional Inputs:
%
%   **subplotrows:**
%        Number of subplot rows (i). Default: 1.
%
%   **subplotcols:**
%        Number of subplot columns (j). Default: 1.
%
%   **keep_existing_axes:**
%        Logical. If true, do NOT clear the figure when reusing or creating
%        it (note: this is the inverse of an internal doclear flag).
%        Default: false (i.e., clear by default).
%
%   **force_resize:**
%        Logical. If true, recompute and apply the size/position even when
%        reusing an existing figure. Default: false. (New figures are
%        always sized and placed.)
%
% :Outputs:
%
%   **f1:**
%        Handle to the figure.
%
%   **axh:**
%        Vector of axes handles for the requested subplots, in subplot
%        index order. Empty when keep_existing_axes is true.
%
%   **tag:**
%        The figure's Tag string (same as the input tagname when supplied).
%
% :Sizing model:
%
%   Per-axis target: ~320 x 260 pixels. Total figure size grows with the
%   number of rows/columns. Both dimensions are clamped to 80% of the
%   primary screen; if the natural size exceeds the cap, the figure is
%   scaled down uniformly so axes stay equal in aspect.
%
% :Placement model:
%
%   For new figures, a grid of candidate top-left positions is scanned
%   (30-pixel step). The first candidate whose bounding box does not
%   overlap any other visible figure is used. If every candidate overlaps,
%   the figure cascades from the top-left by 30 pixels per existing figure.
%
% :Examples:
% ::
%
%    % 1) Single axis (default)
%    [f, ax, tag] = create_figure('one axis');
%    plot(ax, randn(50,1));
%
%    % 2) One row, two columns (1 x 2): side-by-side panels
%    [f, ax] = create_figure('1x2 demo', 1, 2);
%    plot(ax(1), randn(50,1)); title(ax(1), 'left');
%    plot(ax(2), randn(50,1)); title(ax(2), 'right');
%
%    % 3) Two rows, one column (2 x 1): stacked panels
%    [f, ax] = create_figure('2x1 demo', 2, 1);
%    plot(ax(1), randn(50,1)); title(ax(1), 'top');
%    plot(ax(2), randn(50,1)); title(ax(2), 'bottom');
%
%    % 4) One row, three columns (1 x 3): wide strip
%    [f, ax] = create_figure('1x3 demo', 1, 3);
%    for k = 1:3, plot(ax(k), randn(50,1)); title(ax(k), sprintf('col %d',k)); end
%
%    % 5) Three rows, one column (3 x 1): tall strip
%    [f, ax] = create_figure('3x1 demo', 3, 1);
%    for k = 1:3, plot(ax(k), randn(50,1)); title(ax(k), sprintf('row %d',k)); end
%
%    % 6) Three rows, two columns (3 x 2): six panels, larger figure
%    [f, ax] = create_figure('3x2 demo', 3, 2);
%    for k = 1:6, plot(ax(k), randn(50,1)); title(ax(k), sprintf('panel %d',k)); end
%
%    % 7) Reuse the figure but force a resize/replace position
%    [f, ax] = create_figure('3x2 demo', 3, 2, false, true);
%
% :See also:
%   subplot, figure, set, get
%

% ..
%    Programmers' notes:
%    1/2017 Tor Wager - aspect ratio scaled by subplot layout
%    8/2018 Tor Wager - added force-resize-existing flag
%    9/2018 Tor Wager - doresize default off
%    2026   Equal per-axis sizing, 80% screen cap, non-overlapping
%           placement, and tag returned as third output.
% ..

axh = [];

if nargin < 1 || isempty(tagname)
    tagname = 'nmdsfig';
end

doclear = true;
createnew = false;
doresize = false;

if length(varargin) > 2 && ~isempty(varargin{3})
    % keep_existing_axes flag inverts doclear
    doclear = ~varargin{3};
end

if length(varargin) > 3 && ~isempty(varargin{4})
    doresize = varargin{4};
end

old = findobj('Tag', tagname);
old = old( strcmp( get(old, 'Type'), 'figure' ) );

if ~isempty(old)
    % Existing figure with this tag -- reuse

    if length(old) > 1
        close(old(2:end))
        old = old(1);
    end

    if doclear, clf(old); end

    f1 = old;

else
    % Create new

    createnew = true;

    f1 = figure;

    set(f1, 'Tag', tagname, 'Name', tagname, 'color', 'white');
    hold on

end

% activate this figure
figure(f1);

% Parse subplot rows/cols
i = 1;
j = 1;
if ~isempty(varargin)
    if ~isempty(varargin{1}), i = max(1, varargin{1}); end
    if length(varargin) > 1 && ~isempty(varargin{2})
        j = max(1, varargin{2});
    end
end

if doclear

    np = max(1, i * j);

    if np == 1
        axh = gca;
        cla
        set(gca, 'FontSize', 14)
        hold on
    else
        for k = 1:np
            axh(k) = subplot(i, j, k);
            cla;
            set(gca, 'FontSize', 14)
            hold on
        end
        axes(axh(1));
    end

end

if createnew || doresize
    set_figure_position(f1, i, j);
end

tag = get(f1, 'Tag');

end % main function


% =========================================================================
% Sizing
% =========================================================================
function set_figure_position(f1, i, j)
% Size the figure so each subplot is ~constant on-screen, capped at 80%
% of the primary screen; then place it without overlapping existing figs.

screensz = get(0, 'screensize');
screen_w = screensz(3);
screen_h = screensz(4);
max_w    = 0.8 * screen_w;
max_h    = 0.8 * screen_h;

% Per-axis target (medium)
ax_w     = 320;
ax_h     = 260;

% Per-axis "gutter+margin" budget. These approximate the chrome around an
% axes (titles, ticks, labels) so total fig size grows roughly linearly
% with the number of rows/cols.
per_col  = ax_w + 60;
per_row  = ax_h + 70;

% Minimum readable figure size for a 1x1 layout
min_w    = 480;
min_h    = 420;

w = max(min_w, j * per_col + 40);
h = max(min_h, i * per_row + 60);

% Scale down uniformly if either dimension exceeds the 80% cap. Uniform
% scaling preserves the aspect chosen for equal-sized axes.
scale = min(1, min(max_w / w, max_h / h));
w = round(w * scale);
h = round(h * scale);

pos = find_nonoverlapping_position(f1, w, h, screensz);
set(f1, 'Units', 'pixels');
set(f1, 'Position', pos);

end


% =========================================================================
% Placement
% =========================================================================
function pos = find_nonoverlapping_position(f1, w, h, screensz)
% Scan a coarse grid of top-left corners and pick the first slot that
% doesn't intersect any visible figure. Cascade if everything is covered.

screen_w = screensz(3);
screen_h = screensz(4);

others = findobj(0, 'Type', 'figure');
others = others(others ~= f1);

others_pos = zeros(length(others), 4);
for k = 1:length(others)
    % Some figures may use non-pixel Units; normalize.
    u = get(others(k), 'Units');
    set(others(k), 'Units', 'pixels');
    others_pos(k, :) = get(others(k), 'Position');
    set(others(k), 'Units', u);
end

step   = 30;
margin = 30;

% Candidate x ranges left->right, y ranges top->bottom
max_x = max(margin, screen_w - w - margin);
min_y = margin;
max_y = max(margin, screen_h - h - margin);

xs = margin:step:max_x;
ys = max_y:-step:min_y;
if isempty(xs), xs = margin; end
if isempty(ys), ys = max_y;  end

for yy = ys
    for xx = xs
        candidate = [xx yy w h];
        if ~any_overlap(candidate, others_pos)
            pos = candidate;
            return
        end
    end
end

% Fall back: cascade from top-left based on figure count
nfig  = length(others);
xoff  = margin + step * mod(nfig, 10);
yoff  = max(margin, max_y - step * mod(nfig, 10));
pos   = [xoff yoff w h];

end


function tf = any_overlap(a, B)
% a = [x y w h]; B is N x 4. Returns true if a intersects any row of B.
tf = false;
if isempty(B), return; end

ax1 = a(1); ay1 = a(2);
ax2 = a(1) + a(3); ay2 = a(2) + a(4);

for k = 1:size(B, 1)
    bx1 = B(k, 1); by1 = B(k, 2);
    bx2 = bx1 + B(k, 3); by2 = by1 + B(k, 4);
    if ax1 < bx2 && ax2 > bx1 && ay1 < by2 && ay2 > by1
        tf = true;
        return
    end
end

end
