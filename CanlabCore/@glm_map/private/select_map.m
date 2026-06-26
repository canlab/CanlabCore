function [map, which_map, wh_image, varargin] = select_map(obj, varargin)
% Private helper for glm_map: select a result map (and optional image) by keyword.
%
% Pulls an optional leading map-name argument and an optional image index from
% varargin, and returns the corresponding statistic_image result map stored in
% the glm_map object, the resolved map name, the requested image index (or []
% if none given), and the remaining varargin.
%
% which_map (first arg if char and a recognized name) is one of:
%   'betas' | 't' | 'contrast' | 'contrast_t'        (default 't')
%
% The image index may be given either as a bare numeric scalar following the
% map name, or as a 'wh_image'/'image' keyword-value pair. For betas/t this
% indexes regressors; for contrast/contrast_t it indexes contrasts.

which_map = 't';
wh_image  = [];

valid = {'betas', 't', 'contrast', 'contrast_t'};

% Optional leading map name
if ~isempty(varargin) && (ischar(varargin{1}) || isstring(varargin{1})) ...
        && any(strcmpi(char(varargin{1}), valid))
    which_map = lower(char(varargin{1}));
    varargin(1) = [];
end

% Optional image index as a bare numeric scalar
if ~isempty(varargin) && isnumeric(varargin{1}) && isscalar(varargin{1})
    wh_image = varargin{1};
    varargin(1) = [];
end

% Optional image index as a keyword-value pair
wh = find(strcmpi(varargin, 'wh_image') | strcmpi(varargin, 'image'));
if ~isempty(wh)
    wh_image = varargin{wh(1) + 1};
    varargin(wh(1):wh(1) + 1) = [];
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
