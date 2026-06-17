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

[map, which_map, wh_image, varargin] = select_map(obj, varargin{:});

if isempty(map)
    error('glm_map:NoMap', 'Requested map ''%s'' is empty. Fit the model first.', which_map);
end

% statistic_image.table requires a single image; resolve one.
nimg = size(map.dat, 2);
if nimg > 1
    if isempty(wh_image), wh_image = 1; end
    labels = map.image_labels;
    if iscell(labels) && numel(labels) >= wh_image && ~isempty(labels{wh_image})
        thislabel = labels{wh_image};
    else
        thislabel = sprintf('image %d', wh_image);
    end
    fprintf('  %s map has %d images; showing %s (%s). Pass an index to choose another, e.g. table(obj, ''%s'', k).\n', ...
        which_map, nimg, num2str(wh_image), thislabel, which_map);
    map = get_wh_image(map, wh_image);
elseif ~isempty(wh_image)
    map = get_wh_image(map, wh_image);
end

[varargout{1:nargout}] = table(map, varargin{:});

end % table
