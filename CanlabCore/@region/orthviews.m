function orthviews(obj, varargin)
% orthviews Display a region object on SPM orthviews via cluster_orthviews.
%
% Forward the region object and all optional inputs to
% cluster_orthviews and apply the standard CANlab hot/cool colormap to
% the SPM orthviews window.
%
% :Usage:
% ::
%
%     orthviews(obj, [optional inputs])
%
% :Inputs:
%
%   **obj:**
%        A region-class object array.
%
% :Optional Inputs:
%
%   Any inputs accepted by cluster_orthviews.
%
% :Examples:
% ::
%
%     % Match colors left/right
%     all_colors = match_colors_left_right(obj);
%     orthviews(r, all_colors);
%
% :See also:
%   - cluster_orthviews
%   - match_colors_left_right
%   - spm_orthviews_hotcool_colormap

cluster_orthviews(obj, varargin{:});

try
    spm_orthviews_hotcool_colormap([], 0);
catch
end

end