function obj = create_orthogonal_contrast_set(obj, varargin)
% Assign an orthogonal set of contrasts spanning the regressors of interest.
%
% Uses create_orthogonal_contrast_set.m to build a set of mutually orthogonal,
% zero-sum contrasts that span the space of differences among the regressors
% of interest (obj.wh_interest), and stores them in obj.contrasts with
% placeholder names. Nuisance covariates and the intercept get zero weight.
% Any previously assigned contrasts are replaced.
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
%        A glm_map object whose regressors of interest are defined
%        (obj.wh_interest has at least 2 true entries).
%
% :Optional Inputs:
%
%   **'names', {cellstr}:**
%        Names for the contrasts (default 'OrthC1', 'OrthC2', ...).
%
% :Outputs:
%
%   **obj:**
%        glm_map with obj.contrasts ([regressors x contrasts]) and
%        obj.contrast_names set to the orthogonal set.
%
% :Examples:
% ::
%
%     d = fmri_glm_design_matrix(2,'nscan',120,'units','secs', ...
%             'onsets',{[10 40]' [25 55]' [12 42]'}, 'condition_names',{'A','B','C'});
%     g = glm_map(d); g = build_design(g);
%     g = create_orthogonal_contrast_set(g);   % 2 orthogonal contrasts over A,B,C
%     g = diagnostics(g);
%
% :See also:
%   - create_orthogonal_contrast_set (function), add_contrasts, diagnostics, calcEfficiency
%
% ..
%    2026 - Initial implementation.
% ..

names = {};
wh = find(strcmpi(varargin, 'names'));
if ~isempty(wh), names = varargin{wh(1) + 1}; end

whI = find(obj.wh_interest);

if isempty(whI)
    error('glm_map:NoInterestRegressors', '%s', local_interest_guidance());
end
if numel(whI) < 2
    error('glm_map:TooFewInterest', ...
        ['Orthogonal contrasts require at least 2 regressors of interest; found %d. ' ...
         'Add more event types of interest, or check obj.wh_interest.'], numel(whI));
end

% Orthogonal contrast set over the of-interest regressors (rows = contrasts)
Cint = create_orthogonal_contrast_set(numel(whI));   % [(n-1) x n], dispatches to the function

% Embed into the full design width (zeros on nuisance / intercept columns)
ncon = size(Cint, 1);
Cfull = zeros(ncon, obj.num_regressors);
Cfull(:, whI) = Cint;

if isempty(names)
    names = arrayfun(@(i) sprintf('OrthC%d', i), 1:ncon, 'UniformOutput', false);
elseif numel(names) ~= ncon
    error('glm_map:ContrastNames', 'Number of names (%d) must match number of contrasts (%d).', numel(names), ncon);
end

% Store (replace any existing contrasts). glm_map stores [regressors x contrasts].
obj.contrasts      = Cfull';
obj.contrast_names = names(:)';

obj.history{end + 1} = sprintf(['create_orthogonal_contrast_set: %d orthogonal contrast(s) ' ...
    'spanning %d regressors of interest'], ncon, numel(whI));

end % create_orthogonal_contrast_set


% =========================================================================
function s = local_interest_guidance()
s = sprintf([ ...
    'No regressors of interest are defined for this glm_map.\n' ...
    'Regressors of interest are the task-event regressors.\n' ...
    '  - Event / 1st-level designs: add event onsets and build the design\n' ...
    '    (g = build_design(g), or import_onsets / import_SPM); event\n' ...
    '    regressors are then flagged as of interest automatically.\n' ...
    '  - Direct / group designs: by default all non-intercept columns are of\n' ...
    '    interest. Mark covariates of no interest with\n' ...
    '    g.nuisance_columns = [column indices]; at least one column must\n' ...
    '    remain of interest.']);
end
