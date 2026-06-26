function obj = squeeze_figure(obj, varargin)
% Remove top/bottom white space from montage figures, keeping slice proportions.
%
% A montage figure is usually much taller than the slices it contains (the
% slices occupy a horizontal band in the middle, with white space above and
% below). squeeze_figure shrinks each montage figure's height and re-positions
% the slice axes to fill it, scaled so that each slice keeps its on-screen
% pixel size, aspect ratio, and horizontal placement — only the empty vertical
% margins are removed.
%
% Figures that also contain surface views are skipped (resizing would distort
% the surfaces); use this on montage-only figures.
%
% :Usage:
% ::
%
%     o2 = canlab_results_fmridisplay(t);
%     squeeze_figure(o2);                 % tighten the montage figure(s)
%     squeeze_figure(o2, 'margin', 0.05); % leave a 5% top/bottom margin
%
% :Inputs:
%
%   **obj:** an fmridisplay object (handle) with one or more montages.
%
% :Optional Inputs:
%
%   **'margin':** fraction of figure height to leave above and below the
%                 slices (default 0.02).
%
% :Outputs:
%
%   **obj:** the same handle (figures resized in place).
%
% :See also:
%   - montage, canlab_results_fmridisplay, squeeze_axes
%
% ..
%    2026 visualization overhaul
% ..

marg = 0.02;
wh = find(strcmp(varargin, 'margin'), 1);
if ~isempty(wh), marg = varargin{wh + 1}; end
marg = max(0, min(0.4, marg));

% Collect montage slice axes and the figures they live in
mont_ax = gobjects(0);
for m = 1:numel(obj.montage)
    ah = obj.montage{m}.axis_handles;
    mont_ax = [mont_ax, ah(ishandle(ah))]; %#ok<AGROW>
end
if isempty(mont_ax), return, end

mont_figs = arrayfun(@(a) ancestor(a, 'figure'), mont_ax);

% Figures that also contain surfaces — skip those (don't distort surfaces)
surf_figs = gobjects(0);
for s = 1:numel(obj.surface)
    h = obj.surface{s}.object_handle; h = h(ishandle(h));
    if ~isempty(h), surf_figs(end + 1) = ancestor(h(1), 'figure'); end %#ok<AGROW>
end

for f = unique(mont_figs)
    if ~isempty(surf_figs) && any(surf_figs == f), continue, end

    axf = mont_ax(mont_figs == f);
    P = get(axf, 'Position');
    if iscell(P), P = cell2mat(P); end           % N x 4 [x y w h], normalized

    ymin = min(P(:, 2));
    ymax = max(P(:, 2) + P(:, 4));
    cf   = ymax - ymin;                          % content height fraction
    if cf <= 0 || cf >= (1 - 2 * marg), continue, end   % nothing to gain

    scale = (1 - 2 * marg) / cf;

    % Remap each slice axis to fill [marg, 1-marg] vertically (x untouched)
    for a = axf
        p = a.Position;
        a.Position = [p(1), marg + (p(2) - ymin) * scale, p(3), p(4) * scale];
    end

    % Shrink the figure height by the same factor so slices keep their pixel size
    f.Position(4) = round(f.Position(4) * cf / (1 - 2 * marg));
end

if ~iscell(obj.history), obj.history = {}; end
obj.history{end + 1} = 'squeeze_figure: removed vertical white space from montage figure(s)';

end
