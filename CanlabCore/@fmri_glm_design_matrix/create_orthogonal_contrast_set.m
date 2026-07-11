function obj = create_orthogonal_contrast_set(obj, varargin)
% Assign an orthogonal set of contrasts spanning the regressors of interest.
%
% Uses create_orthogonal_contrast_set.m to build a set of mutually orthogonal,
% zero-sum contrasts spanning the space of differences among the regressors of
% interest (the H partition, obj.xX.iH), and stores them in obj.xX.contrasts
% (one row per contrast, one column per design regressor) with placeholder
% names in obj.xX.contrast_names. Nuisance covariates and session constants
% get zero weight. The model must be built first (so obj.xX.iH is defined).
%
% :Usage:
% ::
%
%     obj = create_orthogonal_contrast_set(obj)
%     obj = create_orthogonal_contrast_set(obj, 'names', {'c1','c2',...})
%
% :Inputs:
%
%   **obj:**
%        A built fmri_glm_design_matrix object (obj.xX.iH non-empty).
%
% :Optional Inputs:
%
%   **'names', {cellstr}:**  contrast names (default 'OrthC1', ...).
%
% :Outputs:
%
%   **obj:**
%        with obj.xX.contrasts [n_contrasts x n_regressors] and
%        obj.xX.contrast_names set.
%
% :See also:
%   - create_orthogonal_contrast_set (function), build, glm_map.create_orthogonal_contrast_set
%
% ..
%    2026 - Initial implementation.
% ..

names = {};
wh = find(strcmpi(varargin, 'names'));
if ~isempty(wh), names = varargin{wh(1) + 1}; end

% Regressors of interest = H partition. Requires a built design.
whI = [];
if isstruct(obj.xX) && isscalar(obj.xX) && isfield(obj.xX, 'iH') && ~isempty(obj.xX.X)
    whI = obj.xX.iH(:)';
end

if isempty(whI)
    error('fmri_glm_design_matrix:NoInterestRegressors', '%s', local_interest_guidance());
end
if numel(whI) < 2
    error('fmri_glm_design_matrix:TooFewInterest', ...
        'Orthogonal contrasts require at least 2 regressors of interest (obj.xX.iH); found %d.', numel(whI));
end

nreg = size(obj.xX.X, 2);
Cint = create_orthogonal_contrast_set(numel(whI));   % [(n-1) x n], rows = contrasts
ncon = size(Cint, 1);

C = zeros(ncon, nreg);
C(:, whI) = Cint;

if isempty(names)
    names = arrayfun(@(i) sprintf('OrthC%d', i), 1:ncon, 'UniformOutput', false);
elseif numel(names) ~= ncon
    error('fmri_glm_design_matrix:ContrastNames', 'Number of names (%d) must match number of contrasts (%d).', numel(names), ncon);
end

obj.xX(1).contrasts      = C;
obj.xX(1).contrast_names = names(:)';

if ~iscell(obj.history), obj.history = {}; end
obj.history{end + 1} = sprintf('create_orthogonal_contrast_set: %d contrast(s) over %d interest regressors', ncon, numel(whI));

end % create_orthogonal_contrast_set


% =========================================================================
function s = local_interest_guidance()
s = sprintf([ ...
    'No regressors of interest are defined (obj.xX.iH is empty).\n' ...
    'Regressors of interest are the event/task regressors (the H partition).\n' ...
    'Add event onsets and build the design first, e.g.:\n' ...
    '  d = fmri_glm_design_matrix(TR, ''nscan'', nscan, ''units'', ''secs'', ...\n' ...
    '          ''onsets'', onsets, ''condition_names'', names);\n' ...
    '  d = build(d);\n' ...
    'After build, the event regressors are placed in obj.xX.iH.']);
end
