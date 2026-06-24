function obj = rethreshold(obj, input_threshold, varargin)
% Re-threshold existing fmridisplay blob layers from their retained source.
%
% Recomputes each target layer's displayed voxels from the source object the
% layer retained at addblobs time (see VISUALIZATION_OVERHAUL_NOTES.md), then
% re-renders in place. Because the layer keeps the (un-rethresholded) source,
% you can move the threshold UP or DOWN when the source is a statistic_image /
% image_vector. When the source is only a region, the threshold is an absolute
% cutoff on the stored statistic values (stringency can only increase).
%
% :Usage:
% ::
%
%     obj = rethreshold(obj, p)                       % p-value, 'unc', all layers
%     obj = rethreshold(obj, p, thresh_type)          % thresh_type: 'unc'|'fdr'|'fwe'
%     obj = rethreshold(obj, p, thresh_type, 'k', n)  % cluster-extent k
%     obj = rethreshold(obj, p, ..., 'layers', idx)   % only the listed layers
%
% :Inputs:
%
%   **obj:**
%        An fmridisplay object (handle) with one or more activation_maps.
%
%   **input_threshold:**
%        For statistic_image/image_vector sources, the threshold passed to
%        threshold() (a p-value for 'unc'/'fdr'). For region sources, an
%        absolute magnitude cutoff on the stored statistic.
%
% :Optional Inputs:
%
%   **thresh_type:**
%        'unc' (default), 'fdr', or 'fwe' — only for statistic_image sources.
%
%   **'k':**
%        Followed by a cluster-extent (voxels) — only for statistic_image sources.
%
%   **'layers':**
%        Followed by a numeric vector of layer indices. Default: all layers.
%
% :Outputs:
%
%   **obj:**
%        The same handle, with target layers recomputed and redrawn.
%
% :Examples:
% ::
%
%     t  = threshold(ttest(load_image_set('emotionreg')), .005, 'unc');
%     o2 = canlab_results_fmridisplay(t);
%     o2 = rethreshold(o2, .05,  'unc');   % loosen threshold, blobs grow
%     o2 = rethreshold(o2, .001, 'unc');   % tighten threshold, blobs shrink
%
% :See also:
%   - addblobs, refresh, set_colormap, set_opacity, threshold
%
% ..
%    2026 visualization overhaul
% ..

% Parse options
wh_layers = 1:numel(obj.activation_maps);
whl = find(strcmp(varargin, 'layers'));
if ~isempty(whl)
    wh_layers = varargin{whl(1) + 1};
    varargin(whl(1):whl(1) + 1) = [];
end

thresh_type = 'unc';
if ~isempty(varargin) && (ischar(varargin{1}) || isstring(varargin{1})) ...
        && any(strcmpi(varargin{1}, {'unc', 'fdr', 'fwe'}))
    thresh_type = char(varargin{1});
    varargin(1) = [];
end

k_extent = [];
whk = find(strcmp(varargin, 'k'));
if ~isempty(whk), k_extent = varargin{whk(1) + 1}; end

for k = wh_layers

    if k < 1 || k > numel(obj.activation_maps), continue, end

    layer = obj.activation_maps{k};
    if ~isfield(layer, 'source_object') || isempty(layer.source_object)
        warning('fmridisplay:rethreshold', ...
            'Layer %d retained no source; cannot rethreshold. Re-add via addblobs.', k);
        continue
    end

    src = layer.source_object;

    % Re-derive a region (cl) at the new threshold
    if isa(src, 'statistic_image')

        % P-value based thresholding (uses .p): thresh_type is 'unc'/'fdr'/'fwe'.
        if isempty(k_extent)
            tt = threshold(src, input_threshold, thresh_type, 'noverbose');
        else
            tt = threshold(src, input_threshold, thresh_type, 'k', k_extent, 'noverbose');
        end
        cl = region(tt, 'noverbose');

    elseif isa(src, 'image_vector')

        % Raw-value thresholding (fmri_data / mask / mean image: no p-values).
        % A scalar threshold is treated as a magnitude cutoff (keep |val| > v);
        % a [lo hi] range uses thresh_type ('raw-outside' default, or
        % 'raw-between' if explicitly requested).
        if isscalar(input_threshold)
            rng   = [-abs(input_threshold) abs(input_threshold)];
            rtype = 'raw-outside';
        else
            rng   = input_threshold;
            rtype = thresh_type;
            if ~any(strcmpi(rtype, {'raw-between', 'raw-outside'})), rtype = 'raw-outside'; end
        end
        if isempty(k_extent)
            tt = threshold(src, rng, rtype, 'noverbose');
        else
            tt = threshold(src, rng, rtype, 'k', k_extent, 'noverbose');
        end
        cl = region(tt, 'noverbose');

    elseif isa(src, 'region')

        % Absolute-value cutoff on the stored statistic (raise stringency).
        iv = region2imagevec(src);
        iv.dat(abs(iv.dat) < input_threshold) = 0;
        cl = region(iv, 'noverbose');

    else
        warning('fmridisplay:rethreshold', ...
            'Layer %d source is a %s; rethreshold supports statistic_image/image_vector/region.', ...
            k, class(src));
        continue
    end

    if isempty(cl)
        % Nothing survives: delete this layer's blobs but keep the entry/source.
        if isfield(layer, 'blobhandles') && ~isempty(layer.blobhandles)
            wh = ishandle(layer.blobhandles);
            if any(wh), delete(layer.blobhandles(wh)); end
        end
        obj.activation_maps{k}.blobhandles = [];
        obj.activation_maps{k}.source_region = cl;
        obj.activation_maps{k}.applied_threshold = input_threshold;
        continue
    end

    % Rebuild mapdata + space exactly as addblobs does
    [~, mask] = clusters2mask2011(cl);
    XYZ = cat(2, cl(:).XYZ);
    dim = max(XYZ, [], 2);
    V = struct('mat', cl(1).M, 'dim', dim');
    SPACE = map_to_world_space(V);

    obj.activation_maps{k}.mapdata = mask;
    obj.activation_maps{k}.V = V;
    obj.activation_maps{k}.SPACE = SPACE;
    obj.activation_maps{k}.source_region = cl;
    obj.activation_maps{k}.applied_threshold = input_threshold;
end

obj = refresh(obj, wh_layers);

% Keep an open controller's threshold/colormap/opacity fields in sync.
obj = update_controller(obj);

end
