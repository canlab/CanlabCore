function obj = set_colormap(obj, varargin)
% Change the colors of existing fmridisplay blob layers and re-render in place.
%
% Replaces the color options on already-rendered blob layers without the
% caller re-supplying the source. Accepts the same color keywords as addblobs
% / render_blobs ('color', 'maxcolor'/'mincolor', 'splitcolor', 'colormap',
% 'cmaprange', 'onecolor'/'solid'), strips any conflicting prior color options
% from each target layer, and calls refresh (handle class; see
% VISUALIZATION_OVERHAUL_NOTES.md).
%
% :Usage:
% ::
%
%     obj = set_colormap(obj, 'color', [1 0 0])                 % solid red
%     obj = set_colormap(obj, 'maxcolor', [1 1 0], 'mincolor', [1 0 0])
%     obj = set_colormap(obj, 'splitcolor', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]})
%     obj = set_colormap(obj, 'maxcolor', [1 1 0], 'mincolor', [1 0 0], 'cmaprange', [0 5])
%     obj = set_colormap(obj, ..., 'layers', idx)               % target specific layers
%
% :Inputs:
%
%   **obj:**
%        An fmridisplay object (handle) with one or more activation_maps.
%
% :Optional Inputs:
%
%   **color keywords:**
%        Any color options understood by render_blobs (see addblobs help).
%
%   **'layers':**
%        Followed by a numeric vector of layer indices. Default: all layers.
%
% :Outputs:
%
%   **obj:**
%        The same handle, with target layers re-rendered using the new colors.
%
% :See also:
%   - addblobs, refresh, set_opacity, rethreshold
%
% ..
%    2026 visualization overhaul
% ..

% Pull out 'layers' (our control option) from the color options to splice in
wh_layers = 1:numel(obj.activation_maps);
color_args = varargin;
whl = find(strcmp(color_args, 'layers'));
if ~isempty(whl)
    wh_layers = color_args{whl(1) + 1};
    color_args(whl(1):whl(1) + 1) = [];
end

if isempty(color_args)
    error('fmridisplay:set_colormap:noColor', ...
        'Supply at least one color option, e.g. set_colormap(obj, ''color'', [1 0 0]).');
end

% Color-controlling tokens in render_blobs. Those followed by a value vs flags.
color_keys_val   = {'color', 'maxcolor', 'mincolor', 'splitcolor', 'colormap', 'cmaprange', 'indexmap'};
color_keys_noval = {'onecolor', 'solid', 'indexmap'};

for k = wh_layers

    if k < 1 || k > numel(obj.activation_maps), continue, end
    if ~isfield(obj.activation_maps{k}, 'render_args'), continue, end

    args = obj.activation_maps{k}.render_args;
    if isempty(args), args = {}; end

    % Strip prior color options, then splice in the new ones
    args = strip_keys(args, color_keys_val, true);
    args = strip_keys(args, color_keys_noval, false);
    args = [args, color_args]; %#ok<AGROW>

    obj.activation_maps{k}.render_args = args;
end

obj = refresh(obj, wh_layers);

end


function args = strip_keys(args, keys, has_value)
% Remove keyword tokens (and, if has_value, the following arg) from a cell.
% Flag-style keys (e.g. 'solid') and value-style keys ('color', [..]) of the
% same name are disambiguated by whether the next element is non-char.
i = 1;
out = {};
while i <= numel(args)
    if ischar(args{i}) && any(strcmp(args{i}, keys))
        next_is_value = has_value && (i + 1) <= numel(args) && ~ischar(args{i + 1});
        if next_is_value, i = i + 2; else, i = i + 1; end
        continue
    end
    out{end + 1} = args{i}; %#ok<AGROW>
    i = i + 1;
end
args = out;
end
