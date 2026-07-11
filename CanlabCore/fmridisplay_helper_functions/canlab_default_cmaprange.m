function cmaprange = canlab_default_cmaprange(mapd, varargin)
% Uniform default colour-limit (cmaprange) policy for CanlabCore blob displays.
%
% Single source of truth for the DEFAULT value->colour range used when the user
% does not pass an explicit 'cmaprange'. Montages (render_blobs), volumetric
% surface rendering (render_on_surface, used by fmri_data.surface), and the
% stateful fmridisplay surface path (render_layer_surfaces) all call this, so a
% montage and a surface of the same data scale identically by default.
%
% The policy is robust (percentile-based, NOT [min max]) so that one or a few
% extreme voxels do not compress the colour scale:
%
%   - signed, non-split   : [prctile 2,  prctile 98]
%   - single-sign         : [prctile 10, prctile 90]
%   - 'splitcolor'        : per-arm 4-element range
%                           [p20(neg) p80(neg) p20(pos) p80(pos)], so the
%                           positive and negative colour arms are scaled
%                           independently. The split percentile relaxes
%                           (20 -> 15 -> ...) if the data lack the variability
%                           to give four distinct edges.
%
% (This is the long-standing montage default, factored out of
% render_blobs.get_default_cmaprange so surfaces share it.)
%
% :Usage:
% ::
%
%     cmaprange = canlab_default_cmaprange(datavec)
%     cmaprange = canlab_default_cmaprange(datavec, 'splitcolor')
%
% :Inputs:
%
%   **mapd:**
%        Numeric vector (or array) of map values. Zeros, NaNs and Infs are
%        handled internally (zeros/NaNs dropped; Infs clamped to the max finite
%        magnitude).
%
% :Optional Inputs:
%
%   **'splitcolor':**
%        Return the per-arm 4-element split range instead of a 2-element range.
%
% :Outputs:
%
%   **cmaprange:**
%        2-element [lo hi], or 4-element [negmin negmax posmin posmax] for
%        'splitcolor'. (canlab_colormap and render_on_surface both accept either
%        form; render_on_surface collapses a 4-element range to its outer
%        [lo hi] saturation bounds.)
%
% :See also:
%   - render_blobs, render_on_surface, render_layer_surfaces, canlab_colormap

mapd = double(mapd(:));
mapd = mapd(mapd ~= 0 & ~isnan(mapd));          % drop zeros and NaNs

% Edge case: no usable data
if isempty(mapd)
    cmaprange = [0, 1];
    return
end

% Edge case: Inf values -> clamp to max finite magnitude
if any(isinf(mapd))
    whinf = isinf(mapd);
    mapd(whinf) = sign(mapd(whinf)) .* max(abs(mapd(~whinf)));
end

% Edge case: constant data
if numel(unique(mapd)) == 1
    constant_value = unique(mapd);
    cmaprange = [constant_value - 0.1, constant_value + 0.1];
    return
end

% Default for single-ramp / solid maps
cmaprange = double([prctile(mapd, 10), prctile(mapd, 90)]);

% Signed single-ramp: widen (robustly) through zero so negatives are not all
% clamped to the min colour. Positive-only / negative-only keep [10 90] above.
if ~any(strcmp(varargin, 'splitcolor')) && any(mapd < 0) && any(mapd > 0)
    cmaprange = double([prctile(mapd, 2), prctile(mapd, 98)]);
end

% Split maps: scale the positive and negative arms independently.
prct_splitcolor = 20;
if any(strcmp(varargin, 'splitcolor')) && ~any(strcmp(varargin, 'cmaprange'))

    cmaprange = double([ ...
        safe_prctile(mapd(mapd < 0), prct_splitcolor), ...
        safe_prctile(mapd(mapd < 0), 100 - prct_splitcolor), ...
        safe_prctile(mapd(mapd > 0), prct_splitcolor), ...
        safe_prctile(mapd(mapd > 0), 100 - prct_splitcolor)]);

    % Relax the percentile if four distinct edges are not available
    while numel(unique(cmaprange)) < 4
        prct_splitcolor = prct_splitcolor - 5;
        if prct_splitcolor <= 0
            cmaprange([2, 3]) = cmaprange([1, 4]) * 0.9;   % compress middle
            break
        end
        cmaprange = double([ ...
            safe_prctile(mapd(mapd < 0), prct_splitcolor), ...
            safe_prctile(mapd(mapd < 0), 100 - prct_splitcolor), ...
            safe_prctile(mapd(mapd > 0), prct_splitcolor), ...
            safe_prctile(mapd(mapd > 0), 100 - prct_splitcolor)]);
    end
end

end % main function


function p = safe_prctile(data, percentile)
% Percentile that tolerates empty input (returns 0).
if isempty(data)
    p = 0;
else
    p = prctile(data, percentile);
end
end
