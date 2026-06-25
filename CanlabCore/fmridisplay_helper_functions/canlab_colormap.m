classdef canlab_colormap
% canlab_colormap Canonical value->colour mapping for CANlab displays.
%
% A single source of truth that turns data values into RGB colours for ONE
% colormap type (single ramp / split +-/- / solid / indexed), independent of
% the renderer. Montage slices, surfaces, and legends are meant to all consume
% this same object so they agree exactly (see VISUALIZATION_OVERHAUL_NOTES.md,
% "one central value->colour mapping, two renderers").
%
% The mapping matches render_blobs:
%   - single : interpolate mincolor (at range(1)) -> maxcolor (at range(2)).
%   - split  : positive values interpolate minpos (near 0) -> maxpos (extreme);
%              negative values interpolate maxneg (near 0) -> minneg (extreme).
%   - solid  : one colour for every (in-blob) value.
%   - indexed: round(value) indexes into an n x 3 colormap.
%
% :Usage:
% ::
%
%     cm  = canlab_colormap.single([1 0 0], [1 1 0], [0 5]);     % red->yellow over 0..5
%     cm  = canlab_colormap.split([0 0 1],[0 1 1],[1 .5 0],[1 1 0], [-4 4]);
%     cm  = canlab_colormap.solid([1 0 0]);
%     cm  = canlab_colormap.indexed(hot(64));
%     cm  = canlab_colormap.from_render_args(layer.render_args, layer.cmaprange);
%
%     rgb           = map(cm, values);        % N x 3 RGB for a vector of values
%     [vals, vrgb]  = legend_samples(cm, 100) % one colorbar's worth (spans zero for split)
%     lut           = lut(cm, 256);           % 256 x 3 lookup table over the range
%
% :Properties:
%
%   **type:**   'single' | 'split' | 'solid' | 'indexed'
%   **range:**  value range(s): single/solid -> [lo hi];
%               split -> [negmin negmax posmin posmax]; indexed -> [1 nrows]
%   **colors:** cell of colours (single {min,max}; split {minneg,maxneg,minpos,maxpos};
%               solid {color}) or the n x 3 matrix (indexed)
%
% ..
%    2026 visualization overhaul — central colour pipeline
% ..

    properties
        type   {mustBeMember(type, {'single','split','solid','indexed'})} = 'single'
        range  = [0 1]
        colors = {[0 0 1] [1 1 0]}
    end

    methods
        function obj = canlab_colormap(type, colors, range)
            if nargin == 0, return, end
            obj.type = type;
            obj.colors = colors;
            if nargin >= 3, obj.range = range; end
        end

        function rgb = map(obj, v)
            % MAP  N x 3 RGB for a vector of values (NaN where uncoloured).
            v = double(v(:));
            n = numel(v);
            rgb = nan(n, 3);

            switch obj.type
                case 'solid'
                    rgb = repmat(rowcolor(obj.colors{1}), n, 1);

                case 'single'
                    lo = obj.range(1); hi = obj.range(2);
                    w = clamp01((v - lo) ./ nonzero(hi - lo));
                    mn = rowcolor(obj.colors{1}); mx = rowcolor(obj.colors{2});
                    rgb = (1 - w) .* mn + w .* mx;

                case 'split'
                    r = obj.range;                                   % [negmin negmax posmin posmax]
                    minneg = rowcolor(obj.colors{1}); maxneg = rowcolor(obj.colors{2});
                    minpos = rowcolor(obj.colors{3}); maxpos = rowcolor(obj.colors{4});
                    posrange = nonzero(r(4) - r(3));
                    negspan  = nonzero(r(1) - r(2));
                    pos = v > 0; neg = v < 0;
                    if any(pos)
                        % positive: minpos at r(3) (near 0) -> maxpos at r(4) (extreme)
                        wp = clamp01((max(v(pos), r(3)) - r(3)) ./ posrange);
                        rgb(pos, :) = (1 - wp) .* minpos + wp .* maxpos;
                    end
                    if any(neg)
                        % negative: maxneg at r(2) (near 0) -> minneg at r(1) (extreme)
                        wn = clamp01((min(v(neg), r(2)) - r(2)) ./ negspan);
                        rgb(neg, :) = (1 - wn) .* maxneg + wn .* minneg;
                    end

                case 'indexed'
                    cmap = obj.colors;
                    idx = round(v); idx = min(max(idx, 1), size(cmap, 1));
                    rgb = cmap(idx, :);
            end

            % Clamp to [0 1] but preserve NaN (uncoloured) entries — note MATLAB's
            % max(NaN,0) returns 0, which would otherwise turn uncoloured -> black.
            nanmask = isnan(rgb);
            rgb = min(max(rgb, 0), 1);
            rgb(nanmask) = NaN;
        end

        function [vals, rgb] = legend_samples(obj, n)
            % LEGEND_SAMPLES  n sample values across the full range + their colours,
            % for a SINGLE colorbar (split spans zero in one continuous bar).
            if nargin < 2, n = 100; end
            switch obj.type
                case 'indexed'
                    vals = (1:size(obj.colors, 1))';
                    rgb  = obj.map(vals);
                case {'single', 'solid'}
                    vals = linspace(obj.range(1), obj.range(end), n)';
                    rgb  = obj.map(vals);
                case 'split'
                    vals = linspace(obj.range(1), obj.range(4), n)';   % extreme neg -> extreme pos
                    rgb  = obj.map(vals);
            end
            % grey any uncoloured samples (e.g. exactly zero)
            bad = any(isnan(rgb), 2);
            rgb(bad, :) = repmat([.5 .5 .5], nnz(bad), 1);
        end

        function out = lut(obj, n)
            % LUT  n x 3 lookup table spanning the full range (for indexed / true-colour use).
            if nargin < 2, n = 256; end
            [~, out] = obj.legend_samples(n);
        end

        function rgb = colorbar_ramp(obj, n)
            % COLORBAR_RAMP  n x 3 CONTINUOUS colour ramp for a legend bar / preview
            % swatch (low -> high), with NO threshold gap — unlike map()/legend_samples,
            % which grey the sub-threshold middle of a split map. For split, the first
            % half is the negative ramp (minneg -> maxneg) and the second the positive
            % ramp (minpos -> maxpos), so the bar reads extreme-neg ... 0 ... extreme-pos.
            % This is the single source for the controller stripe and the figure legend.
            if nargin < 2, n = 64; end
            switch obj.type
                case 'solid'
                    rgb = repmat(rowcolor(obj.colors{1}), n, 1);
                case 'single'
                    rgb = ramp_between(obj.colors{1}, obj.colors{2}, n);
                case 'split'
                    h = floor(n / 2);
                    rgb = [ramp_between(obj.colors{1}, obj.colors{2}, h); ...
                           ramp_between(obj.colors{3}, obj.colors{4}, n - h)];
                case 'indexed'
                    rgb = obj.colors;
            end
            rgb = min(max(rgb, 0), 1);
        end

        function tf = isequal_to(obj, other)
            tf = isa(other, 'canlab_colormap') && strcmp(obj.type, other.type) && ...
                isequal(obj.colors, other.colors) && isequal(obj.range, other.range);
        end
    end

    methods (Static)
        function obj = single(mincolor, maxcolor, range)
            if nargin < 3 || isempty(range), range = [0 1]; end
            obj = canlab_colormap('single', {mincolor, maxcolor}, range(:)');
        end

        function obj = split(minneg, maxneg, minpos, maxpos, range)
            if nargin < 5 || isempty(range), range = [-1 1]; end
            range = expand_split_range(range);
            obj = canlab_colormap('split', {minneg, maxneg, minpos, maxpos}, range);
        end

        function obj = solid(color)
            obj = canlab_colormap('solid', {color}, [0 1]);
        end

        function obj = indexed(cmap)
            obj = canlab_colormap('indexed', cmap, [1 size(cmap, 1)]);
        end

        function obj = from_render_args(args, clim)
            % Build from the render_args a layer stores (splitcolor / maxcolor+
            % mincolor / color / colormap), plus its value range (cmaprange).
            if nargin < 2, clim = []; end
            if isempty(args), args = {}; end
            hask = @(k) any(strcmp(args, k));
            valk = @(k) args{find(strcmp(args, k), 1) + 1};

            if hask('splitcolor')
                sc = valk('splitcolor');                  % {minneg maxneg minpos maxpos}
                obj = canlab_colormap.split(sc{1}, sc{2}, sc{3}, sc{4}, default_clim(clim, [-1 1]));
            elseif hask('color')
                obj = canlab_colormap.solid(valk('color'));
            elseif hask('maxcolor') || hask('mincolor')
                mx = [1 1 0]; mn = [1 0 0];
                if hask('maxcolor'), mx = valk('maxcolor'); end
                if hask('mincolor'), mn = valk('mincolor'); end
                obj = canlab_colormap.single(mn, mx, default_clim(clim, [0 1]));
            elseif hask('colormap') && isnumeric(valk('colormap'))
                obj = canlab_colormap.indexed(valk('colormap'));
            else
                % addblobs default split (mango)
                obj = canlab_colormap.split([.5 0 1], [0 .8 .3], [1 .2 1], [1 1 .3], default_clim(clim, [-1 1]));
            end
        end
    end
end


% ---- local helpers ------------------------------------------------------

function c = rowcolor(c)
c = double(c(:)');
if numel(c) ~= 3, error('canlab_colormap:badColor', 'Colours must be 1x3 RGB.'); end
end

function rgb = ramp_between(c1, c2, n)
% n x 3 linear interpolation from colour c1 to colour c2.
c1 = rowcolor(c1); c2 = rowcolor(c2);
w = linspace(0, 1, n)';
rgb = (1 - w) .* c1 + w .* c2;
end

function w = clamp01(w)
w = min(max(w, 0), 1);
end

function d = nonzero(d)
if d == 0, d = eps; end
end

function r = expand_split_range(range)
% Accept [lo hi] (expand to [lo 0 0 hi]) or a full 4-element [negmin negmax posmin posmax].
range = range(:)';
if numel(range) == 2
    r = [min(range) 0 0 max(range)];
elseif numel(range) == 4
    r = range;
else
    error('canlab_colormap:badRange', 'split range must be 2 or 4 elements.');
end
end

function out = default_clim(clim, fallback)
if isempty(clim), out = fallback; else, out = clim; end
end
