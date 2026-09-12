function obj = legend(obj, varargin)
% legend Create a colorbar legend for an fmridisplay object.
%
% Creates legend(s) for an fmridisplay object and adds legend axis
% handles to obj.activation_maps{:}.legendhandle.
%
% Continuous (value-mapped) layers get a colour ramp with numeric ticks.
% Indexed layers -- atlas / per-region maps added with 'indexmap', e.g. by
% atlas.montage or region.montage -- get a DISCRETE legend: one colour block
% per region index, labelled with the layer's 'labels' when they were given
% (else with the indices). Both kinds are inferred from each layer's stored
% render options, so legend(obj) with no further arguments draws the right
% legend for every layer; the display controller's 'Toggle legend' button
% relies on this. remove_legend(obj) deletes what legend() drew.
%
% Notes: scaleanchors is min and max values for pos and neg range (for
% splitmap).
%
% :Usage:
% ::
%
%     obj = legend(obj, varargin)
%     obj = legend(obj, 'figure')  % new figure
%
% :Inputs:
%
%   **obj:**
%        An fmridisplay object with one or more activation_maps
%        attached.
%
% :Optional Inputs:
%
%   **{'figure', 'newfig'}:**
%        Plot legends in a new figure with larger panels and font.
%        Default: plot small legends in the current figure.
%
%   **'noverbose':**
%        Suppress informational messages (e.g., 'No variability...').
%
%   **'indexmap':**
%        Followed by a colormap (n x 3) for an indexed/parcellation
%        legend (one colour block per row). Overrides the colormap stored
%        on the layer(s). Normally not needed: layers added with
%        'indexmap' are detected automatically.
%
%   **'labels':**
%        Followed by a cell array of label strings, one per colormap
%        row. Only used for indexed legends (a warning is issued
%        otherwise). Overrides labels stored on the layer(s).
%
% :Outputs:
%
%   **obj:**
%        The fmridisplay object with .activation_maps{c}.legendhandle
%        set to the legend axis handle for each activation map.
%
% :Examples:
% ::
%
%     o2 = legend(o2);            % small legends in current figure
%     o2 = legend(o2, 'figure');  % large legends in a new figure
%
%     % Discrete, labelled legend for an atlas layer
%     atl = select_atlas_subset(load_atlas('canlab2024'), {'Thal'});
%     o2  = montage(atl, 'labels', atl.labels);   % stores 'indexmap' + 'labels'
%     o2  = legend(o2);                           % one colour block per region
%     o2  = remove_legend(o2);                    % and take it off again
%
% :See also:
%   - fmridisplay
%   - addblobs
%   - remove_legend
%
% ..
%    Tor Wager
%    8/17/2016 - pkragel updated to accommodate split colormap
%
%    Michael Sun
%    07/29/2024 - Updated to allow for indexmap labelling
%
%    2026 - 'indexmap' / 'labels' are read from each layer's stored
%    render_args, drawn per layer at the standard legend positions, and
%    tracked in legendhandle so remove_legend / the controller's Toggle
%    legend can find and delete them.
% ..

doverbose = true;
donewfig = false;
user_indexmap = [];      % explicit 'indexmap' colormap (overrides layer-stored)
user_labels = {};        % explicit 'labels' (overrides layer-stored)

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'figure' 'newfig'}, donewfig = true;

            case 'noverbose', doverbose = false;

            case 'indexmap'
                user_indexmap = varargin{i+1};

            case 'labels'
                user_labels = varargin{i+1};

            % otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if ~isempty(user_labels) && isempty(user_indexmap) && ~any_indexed_layer(obj)
    warning('''labels'' doesn''t do anything without an ''indexmap'' argument (or an indexmap layer).');
end

if donewfig
    create_figure('legend'); axis off

    mypositions = {[0.100    0.1200    0.80    0.04] ...
        [0.100    0.2800    0.80    0.04] ...
        [0.100    0.46    0.80    0.04] ...
        [0.100    0.64    0.80    0.04] ...
        };

    myfontsize = 18;

else
    mypositions = {[0.05    0.1200    0.20    0.04] ...
        [0.350    0.1200    0.20    0.04] ...
        [0.650    0.1200    0.20    0.04] ...
        [0.05    0.2800    0.20    0.04] ...
        };
    myfontsize = 14;

end

% The legend axes go in the current figure. Callers that want them on a specific
% figure (e.g. the controller drawing onto the montage figure) make it current
% first; new axes are parented explicitly so a later gcf change can't move them.
legfig = gcf;

for c = 1:length(obj.activation_maps)

    currentmap = obj.activation_maps{c};

    if c > length(mypositions)
        disp('Maximum number of legends exceeded. Not plotting remaining legends.');
        break
    end

    % Indexed (atlas / per-region) layer -> discrete legend
    % -------------------------------------------------------------------
    [cmap, labels] = indexed_spec_for_layer(currentmap, user_indexmap, user_labels);

    if ~isempty(cmap)
        pos = mypositions{c};
        if numel(obj.activation_maps) == 1 && size(cmap, 1) > 12
            pos(3) = 0.90;        % a lone atlas layer with many regions: use the full width
        end
        h = draw_indexed_legend(legfig, pos, cmap, labels, myfontsize);
        obj = register_montage_legend(obj, c, h);
        continue
    end

    % Continuous (value-mapped) layer -> colour ramp
    % -------------------------------------------------------------------
    scaleanchors = currentmap.cmaprange; %for text labels

    if isempty(scaleanchors), continue, end

    % adjuts in case there are no values on one end
    scaleanchors(isnan(scaleanchors)) = 0;

    if any(isinf(scaleanchors))
        warning('Some scale anchor values in fmridisplay obj.activation_maps.cmaprange are Inf. Expect erratic behavior/errors.');
    end

    if ~diff(scaleanchors)
        if doverbose, disp('No variability in mapped values. Not plotting legend.'); end
        continue
    end

    h = axes('Parent', legfig, 'Position', mypositions{c});
    obj = register_montage_legend(obj, c, h);
    hold on;

    % fix for multiple blobs
    if size(currentmap.mincolor, 2) > 3 || size(currentmap.maxcolor, 2) > 3
        fprintf('Warning! Extra colors in map...legend will not display correctly with multiple blobs');

        currentmap.mincolor = currentmap.mincolor(:, 1:3);
        currentmap.maxcolor = currentmap.maxcolor(:, 1:3);
    end

    nsteps = 100;

    wvals = linspace(0, 1, nsteps);

    % separate plot for split colormap
    issplitmap = size(currentmap.mincolor, 1) - 1; % zero for one row, 1 for 2+

    % Guard against inconsistent colour fields (e.g. a stale 2-row split mincolor
    % left after switching to a single colormap): a split legend needs 4 colormap
    % anchors (cmaprange). If they aren't present, treat it as a single map rather
    % than indexing past the end of scaleanchors.
    if issplitmap && numel(scaleanchors) < 4
        issplitmap = 0;
    end

    if issplitmap

        [legvals, rgb] = get_split_colormap_values_and_colors(nsteps, scaleanchors, currentmap);

    else
        legvals = linspace(scaleanchors(1), scaleanchors(2), nsteps);

    end

    for i = 2:nsteps

        if ~issplitmap

            fcolor = (1 - wvals(i)) * currentmap.mincolor(1, :) + wvals(i) * currentmap.maxcolor(1, :);

        else

            fcolor = rgb(i, :);

        end

        fcolor(fcolor > 1) = 1;
        fcolor(fcolor < 0) = 0;

        try
            fill([legvals(i-1) legvals(i-1) legvals(i) legvals(i)], [0 1 1 0], fcolor, 'EdgeColor', 'none');

        catch ME
            % (was a `keyboard` debug stop, which hangs unattended / CI runs)
            warning('fmridisplay:legend:fill', 'Problem with legend fill (step %d): %s', i, ME.message);
            break
        end

    end

    if ~issplitmap

        ntick=5;
        myxtick = linspace(scaleanchors(1), scaleanchors(2), ntick);

    else

        ntick=4;
        myxtick = scaleanchors;  % linspace(scaleanchors(1), scaleanchors(4), ntick);

    end

    % kludgy fix for sig digits
    ticklabels = cell(1, ntick);
    for i = 1:ntick
        ticklabels{i} = sprintf('%3.2f', myxtick(i));
    end

    if sum(strcmp(ticklabels, ticklabels{1})) > 1
        for i = 1:ntick
            ticklabels{i} = sprintf('%3.3f', myxtick(i));
        end
    end

    % Relabel if we need more significant digits; we have duplicates
    if sum(strcmp(ticklabels, ticklabels{1})) > 1
        for i = 1:ntick
            ticklabels{i} = sprintf('%3.4f', myxtick(i));
        end
    end

    % remove duplicates and make sure sorted ascending. needed if missing
    % pos or neg values
    [myxtick, wh] = unique(myxtick);
    ticklabels = ticklabels(wh);

    set(gca, 'YTickLabel', '', 'YColor', 'w', 'FontSize', myfontsize);
    axis tight
    set(gca, 'XLim', [min(scaleanchors) max(scaleanchors)], 'XTick', myxtick, 'XTickLabel', ticklabels);

end

end % function



% -------------------------------------------------------------------------
% Subfunctions
% -------------------------------------------------------------------------

function obj = register_montage_legend(obj, c, h)
% Track a montage legend axes on layer c. A previous montage legend for the layer
% (montage_legendhandle) is deleted -- legend() redraws it -- but surface colorbars
% tracked in legendhandle by render_layer_surfaces are kept, so legendhandle holds
% the layer's complete legend set for remove_legend / removeblobs.
lay = obj.activation_maps{c};
keep = gobjects(0);
if isfield(lay, 'legendhandle') && ~isempty(lay.legendhandle)
    keep = lay.legendhandle(ishandle(lay.legendhandle));
end
if isfield(lay, 'montage_legendhandle') && ~isempty(lay.montage_legendhandle)
    old = lay.montage_legendhandle(ishandle(lay.montage_legendhandle));
    keep = keep(~arrayfun(@(x) any(x == old), keep));   % (ismember is not defined for graphics handles)
    delete(old);
end
obj.activation_maps{c}.montage_legendhandle = h;
obj.activation_maps{c}.legendhandle = [reshape(keep, 1, []), h];
end


function tf = any_indexed_layer(obj)
% True if any layer carries an 'indexmap' colormap in its stored render options.
tf = false;
for c = 1:numel(obj.activation_maps)
    if ~isempty(indexed_spec_for_layer(obj.activation_maps{c}, [], {}))
        tf = true;
        return
    end
end
end


function [cmap, labels] = indexed_spec_for_layer(layer, user_indexmap, user_labels)
% Colormap + labels for a discrete legend: explicit inputs first, else the
% 'indexmap' / 'labels' kept in the layer's render options (addblobs stores the
% final option set in render_args). cmap = [] means "not an indexed layer".
cmap = user_indexmap;
labels = user_labels;

args = {};
if isfield(layer, 'render_args') && iscell(layer.render_args), args = layer.render_args; end

if isempty(cmap)
    wh = find(strcmp(args, 'indexmap'), 1);
    if ~isempty(wh) && wh < numel(args) && isnumeric(args{wh + 1}) && size(args{wh + 1}, 2) == 3
        cmap = args{wh + 1};
    end
end

if isempty(cmap)
    labels = {};
    return
end

if isempty(labels)
    wh = find(strcmp(args, 'labels'), 1);
    if ~isempty(wh) && wh < numel(args), labels = args{wh + 1}; end
end

if ~isempty(labels)
    labels = cellstr(labels);   % cell / string array / char matrix -> cell of char
    labels = labels(:)';
end
end


function h = draw_indexed_legend(fig, pos, cmap, labels, fontsize)
% One filled colour block per colormap row, ticks centred in each block. Tick
% labels are the caller's labels (rotated when there are many) or the region
% indices (thinned when there are many). Returns the legend axes handle, which
% the caller tracks in legendhandle so remove_legend / Toggle legend can delete it.
n = size(cmap, 1);

h = axes('Parent', fig, 'Position', pos);
hold(h, 'on');

for i = 1:n
    fill(h, [i-1 i-1 i i], [0 1 1 0], cmap(i, :), 'EdgeColor', 'none');
end

mid = (1:n) - 0.5;          % block centres
rot = 0;
fs  = fontsize;

if ~isempty(labels)
    nlab = min(numel(labels), n);
    if numel(labels) ~= n
        warning('fmridisplay:legend:labelcount', ...
            'Indexed legend: %d labels for %d colormap rows; labelling the first %d.', numel(labels), n, nlab);
    end
    ticks = mid(1:nlab);
    ticklabels = labels(1:nlab);
    if n > 8
        rot = 90;
        fs  = max(6, min(fontsize, round(fontsize * 12 / n)));
    end
else
    keep = 1:n;
    if n > 20, keep = unique(round(linspace(1, n, 10))); end     % thin index ticks
    ticks = mid(keep);
    ticklabels = arrayfun(@num2str, keep, 'UniformOutput', false);
end

set(h, 'XLim', [0 n], 'YLim', [0 1], 'XTick', ticks, 'XTickLabel', ticklabels, ...
    'XTickLabelRotation', rot, 'YTick', [], 'YColor', 'w', 'FontSize', fs, ...
    'Box', 'off', 'TickLength', [0 0]);
end


function [legvals, rgb] = get_split_colormap_values_and_colors(nsteps, scaleanchors, currentmap)
% Values corresponding to colors
legvals = linspace(scaleanchors(1), scaleanchors(2), nsteps./2)';
legvals = [legvals; 0; 0; scaleanchors(3); linspace(scaleanchors(3), scaleanchors(4), nsteps./2)'];

% Colors
rgb1 = zeros(nsteps./2, 3);
rgb2 = zeros(nsteps./2, 3);

for i = 1:3
    rgb1(:, i) = linspace(currentmap.mincolor(2, i), currentmap.maxcolor(2, i), nsteps./2)';
end

rgb1 = [rgb1; ones(3, 3)]; % pad with 1 1 1 white

for i = 1:3
    rgb2(:, i) = linspace(currentmap.mincolor(1, i), currentmap.maxcolor(1, i), nsteps./2)';
end

rgb = [rgb1; rgb2];

end
