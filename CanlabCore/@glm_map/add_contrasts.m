function obj = add_contrasts(obj, C, names, varargin)
% Add one or more linear contrasts to a glm_map object.
%
% Appends rows of contrast weights (over regressors) to obj.contrasts and
% names to obj.contrast_names, validating that the weight vectors match the
% number of regressors in the design.
%
% :Usage:
% ::
%
%     obj = add_contrasts(obj, C, names)
%
% :Inputs:
%
%   **obj:**
%        A glm_map object with a design available (obj.X non-empty).
%
%   **C:**
%        A [num_contrasts x num_regressors] matrix; each ROW is one contrast
%        over the regressors. (Stored internally as [regressors x contrasts].)
%
%   **names:**
%        Cell array of contrast names, one per row of C. Optional; defaults
%        to 'Con1', 'Con2', ...
%
% :Outputs:
%
%   **obj:**
%        glm_map with contrasts and contrast_names appended.
%
% :Examples:
% ::
%
%     g = glm_map('X', [ones(30,1) zscore((1:30)') randn(30,1)], 'level', 2);
%     g = add_contrasts(g, [0 1 0; 0 0 1], {'slope' 'nuisance'});
%
% :See also:
%   - diagnostics, fmri_data.regress
%
% ..
%    Programmers' notes:
%    2026 - Initial implementation.
% ..

if nargin < 3 || isempty(names)
    names = {};
end
if ~iscell(names), names = cellstr(names); end

nreg = obj.num_regressors;
if nreg == 0
    error('glm_map:NoDesign', 'No design available; set obj.X or build a design before adding contrasts.');
end

% Each row of C is a contrast. A contrast may omit the trailing
% intercept/baseline weight(s): if it is one or more elements short, pad with
% 0 weights. (fit also pads contrasts to the final design width, so the
% intercept need not be weighted explicitly.)
if size(C, 2) > nreg
    error('glm_map:ContrastSize', ...
        'Each contrast may have at most %d weights (one per regressor); got %d columns in C.', nreg, size(C, 2));
elseif size(C, 2) < nreg
    C = [C, zeros(size(C, 1), nreg - size(C, 2))];
end

ncon_new = size(C, 1);

% Default names
if isempty(names)
    start = numel(obj.contrast_names);
    names = arrayfun(@(k) sprintf('Con%d', start + k), 1:ncon_new, 'UniformOutput', false);
elseif numel(names) ~= ncon_new
    error('glm_map:ContrastNames', 'Number of names (%d) must match number of contrasts (%d).', numel(names), ncon_new);
end

% Append. Internally contrasts are stored [regressors x contrasts].
obj.contrasts = [obj.contrasts, C'];
obj.contrast_names = [obj.contrast_names(:); names(:)]';

obj.history{end + 1} = sprintf('add_contrasts: added %d contrast(s)', ncon_new);

end % add_contrasts
