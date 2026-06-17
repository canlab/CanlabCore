function varargout = table(obj, varargin)
% Print a table of significant regions for a chosen glm_map result map.
%
% Selects one of the fitted statistic_image maps and delegates to
% statistic_image.table (which forms regions and labels them against an atlas).
%
% :Usage:
% ::
%
%     table(obj)                  % default: thresholded t map
%     table(obj, which_map, ...)  % which_map selects the map
%     [results_table, ...] = table(obj, ...)
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
%        Remaining arguments pass through to statistic_image.table.
%
% :Outputs:
%
%        Whatever statistic_image.table returns (e.g. a results table object).
%
% :Examples:
% ::
%
%     table(g, 'contrast');
%
% :See also:
%   - statistic_image.table, region, montage
%
% ..
%    Programmers' notes:
%    2026 - Initial implementation (thin delegation).
% ..

[map, which_map, varargin] = select_map(obj, varargin{:});

if isempty(map)
    error('glm_map:NoMap', 'Requested map ''%s'' is empty. Fit the model first.', which_map);
end

[varargout{1:nargout}] = table(map, varargin{:});

end % table
