function obj = set_opacity(obj, value, varargin)
% Set blob opacity on existing fmridisplay layers and re-render in place.
%
% Adjusts the transparency of already-rendered blob layers without the caller
% re-supplying the source. Mutates each target layer's stored render options
% and calls refresh (handle class; see VISUALIZATION_OVERHAUL_NOTES.md).
%
% :Usage:
% ::
%
%     obj = set_opacity(obj, value)                  % constant opacity, all layers
%     obj = set_opacity(obj, value, 'layers', idx)   % only the listed layers
%     obj = set_opacity(obj, 'scaled')               % opacity scales with voxel value
%     obj = set_opacity(obj, 'scaled', 'layers', idx)
%
% :Inputs:
%
%   **obj:**
%        An fmridisplay object (handle) with one or more activation_maps.
%
%   **value:**
%        Opacity in [0 1] (1 = opaque, 0 = invisible), OR the string
%        'scaled' to map opacity to voxel magnitude.
%
% :Optional Inputs:
%
%   **'layers':**
%        Followed by a numeric vector of layer indices. Default: all layers.
%
% :Outputs:
%
%   **obj:**
%        The same handle, with target layers re-rendered at the new opacity.
%
% :Examples:
% ::
%
%     o2 = canlab_results_fmridisplay(t);
%     o2 = set_opacity(o2, 0.4);          % make all blobs 40% opaque
%     o2 = set_opacity(o2, 'scaled');     % opacity follows voxel value
%
% :See also:
%   - addblobs, refresh, set_colormap, rethreshold
%
% ..
%    2026 visualization overhaul
% ..

wh_layers = 1:numel(obj.activation_maps);
whl = find(strcmp(varargin, 'layers'));
if ~isempty(whl), wh_layers = varargin{whl(1) + 1}; end

% Tokens that control transparency in render_blobs (and their trailing values)
trans_keys_noval = {'trans', 'transparent', 'scaledtransparency'};
trans_keys_val   = {'constanttrans', 'transvalue'};

for k = wh_layers

    if k < 1 || k > numel(obj.activation_maps), continue, end
    if ~isfield(obj.activation_maps{k}, 'render_args'), continue, end

    args = obj.activation_maps{k}.render_args;
    if isempty(args), args = {}; end

    % Strip any existing transparency options
    args = strip_keys(args, trans_keys_noval, false);
    args = strip_keys(args, trans_keys_val, true);

    % Append the requested transparency
    if (ischar(value) || isstring(value)) && strcmpi(value, 'scaled')
        args{end + 1} = 'scaledtransparency'; %#ok<AGROW>
    else
        validateattributes(value, {'numeric'}, {'scalar', '>=', 0, '<=', 1}, mfilename, 'value');
        args{end + 1} = 'transvalue'; %#ok<AGROW>
        args{end + 1} = double(value); %#ok<AGROW>
    end

    obj.activation_maps{k}.render_args = args;
end

obj = refresh(obj, wh_layers);
obj = update_controller(obj);

end


function args = strip_keys(args, keys, has_value)
% Remove keyword tokens (and, if has_value, the following arg) from a cell.
i = 1;
out = {};
while i <= numel(args)
    if ischar(args{i}) && any(strcmp(args{i}, keys))
        if has_value, i = i + 2; else, i = i + 1; end
        continue
    end
    out{end + 1} = args{i}; %#ok<AGROW>
    i = i + 1;
end
args = out;
end
