function orthviews(obj, varargin)
% orthviews(obj, varargin)
%
% Uses cluster_orthviews.m
% Takes all inputs to cluster_orthviews.m
%
% Examples:
%
% - match colors left/right
% all_colors = match_colors_left_right(obj);
% orthviews(r, all_colors);

cluster_orthviews(obj, varargin{:});

end