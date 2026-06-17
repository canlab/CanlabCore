function o2 = montage(obj, varargin)
% Display a montage of a chosen result map from a glm_map object.
%
% Selects one of the fitted statistic_image maps (betas, t, contrast
% estimates, or contrast t) and delegates to statistic_image/image_vector
% montage for display.
%
% :Usage:
% ::
%
%     montage(obj)                 % default: thresholded t map
%     montage(obj, which_map, ...) % which_map selects the map
%     o2 = montage(obj, ...)       % return the fmridisplay handle
%
% :Inputs:
%
%   **obj:**
%        A fitted glm_map object.
%
% :Optional Inputs:
%
%   **which_map:**
%        One of 'betas' | 't' | 'contrast' | 'contrast_t' (default 't').
%        Any remaining arguments are passed through to statistic_image.montage.
%
% :Outputs:
%
%   **o2:**
%        The fmridisplay object created by the underlying montage call.
%
% :Examples:
% ::
%
%     montage(g, 't');
%     montage(g, 'contrast', 'trans_white');
%
% :See also:
%   - statistic_image, image_vector.montage, table, plot_design
%
% ..
%    Programmers' notes:
%    2026 - Initial implementation (thin delegation).
% ..

[map, which_map, varargin] = select_map(obj, varargin{:});

if isempty(map)
    error('glm_map:NoMap', 'Requested map ''%s'' is empty. Fit the model first.', which_map);
end

if nargout > 0
    o2 = montage(map, varargin{:});
else
    montage(map, varargin{:});
end

end % montage
