function obj = check_properties(obj, varargin)
% Validate and enforce types on a glm_map object's properties.
%
% Lightweight consistency checks: ensures cell-array metadata fields are
% cells, the level is 1 or 2, contrast bookkeeping is consistent, and any
% contained statistic_image maps pass their own check_properties.
%
% :Usage:
% ::
%
%     obj = check_properties(obj)
%
% :Inputs:
%
%   **obj:**
%        A glm_map object.
%
% :Outputs:
%
%   **obj:**
%        The glm_map object with types/fields coerced where needed. Warns on
%        inconsistencies it cannot silently fix.
%
% :See also:
%   - statistic_image.check_properties
%
% ..
%    Programmers' notes:
%    2026 - Initial implementation.
% ..

% Cell-array fields
if ~iscell(obj.contrast_names), obj.contrast_names = cellstr(obj.contrast_names); end
if ~iscell(obj.warnings),       obj.warnings = {};  end
if ~iscell(obj.history),        obj.history  = {};  end
if ~iscell(obj.regressor_names_direct) && ~isempty(obj.regressor_names_direct)
    obj.regressor_names_direct = cellstr(obj.regressor_names_direct);
end

% Level
if ~ismember(obj.level, [1 2])
    warning('glm_map:BadLevel', 'level should be 1 or 2; got %s.', num2str(obj.level));
end

% is_timeseries should be logical
obj.is_timeseries = logical(obj.is_timeseries);

% AR only meaningful for timeseries
if obj.level == 2 && obj.is_timeseries
    warning('glm_map:LevelTimeseries', 'is_timeseries is true but level is 2 (group); AR models are not appropriate for group images.');
end

% Contrast bookkeeping
if ~isempty(obj.contrasts)
    if size(obj.contrasts, 1) ~= obj.num_regressors && obj.num_regressors > 0
        warning('glm_map:ContrastSize', 'contrasts has %d rows but design has %d regressors.', ...
            size(obj.contrasts, 1), obj.num_regressors);
    end
    if ~isempty(obj.contrast_names) && numel(obj.contrast_names) ~= size(obj.contrasts, 2)
        warning('glm_map:ContrastNames', 'Number of contrast_names (%d) does not match number of contrasts (%d).', ...
            numel(obj.contrast_names), size(obj.contrasts, 2));
    end
end

% Delegate to contained statistic_image maps where possible
for f = {'betas', 't', 'contrast_estimates', 'contrast_t'}
    m = obj.(f{1});
    if ~isempty(m) && isa(m, 'statistic_image') && ismethod(m, 'check_properties')
        obj.(f{1}) = check_properties(m);
    end
end

end % check_properties
