function [map, which_map, varargin] = select_map(obj, varargin)
% Private helper for glm_map: select a result map by keyword.
%
% Pulls an optional leading map-name argument from varargin and returns the
% corresponding statistic_image (or fmri_data) result map stored in the
% glm_map object, the resolved name, and the remaining varargin.
%
% which_map (first arg if char and a recognized name) is one of:
%   'betas' | 't' | 'contrast' | 'contrast_t'
% Default is 't'.

which_map = 't';

valid = {'betas', 't', 'contrast', 'contrast_t'};

if ~isempty(varargin) && (ischar(varargin{1}) || isstring(varargin{1})) ...
        && any(strcmpi(char(varargin{1}), valid))
    which_map = lower(char(varargin{1}));
    varargin(1) = [];
end

switch which_map
    case 'betas'
        map = obj.betas;
    case 't'
        map = obj.t;
    case 'contrast'
        map = obj.contrast_estimates;
    case 'contrast_t'
        map = obj.contrast_t;
    otherwise
        error('glm_map:BadMap', 'Unknown map ''%s''. Use betas | t | contrast | contrast_t.', which_map);
end

end % select_map
