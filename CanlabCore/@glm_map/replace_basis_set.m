function obj = replace_basis_set(obj, condition_num, xBF_hires, varargin)
% Replace the basis set for one condition of an event/1st-level glm_map.
%
% Delegates to fmri_glm_design_matrix.replace_basis_set, rebuilds the design
% matrix, and clears any fitted results / contrasts / diagnostics that were
% computed with the previous basis set (they no longer apply once the number
% and meaning of the design columns change). If a data object is supplied (or
% has been retained), the model is re-fit; otherwise design diagnostics are
% recomputed on the new design.
%
% :Usage:
% ::
%
%     obj = replace_basis_set(obj, condition_num, xBF_hires)
%     obj = replace_basis_set(obj, condition_num, xBF_hires, 'data', fmri_data_obj, ...)
%
% :Inputs:
%
%   **obj:**
%        A glm_map object built around an fmri_glm_design_matrix (event /
%        1st-level mode; obj.design non-empty).
%
%   **condition_num:**
%        Index of the condition whose basis set to replace.
%
%   **xBF_hires:**
%        A basis-set struct (e.g. from fmri_spline_basis or spm_get_bf) with a
%        .bf matrix and .dt sampling interval (seconds/sample).
%
% :Optional Inputs:
%
%   **'data', fmri_data_obj:**
%        If provided, re-fit the model on this data after rebuilding. Any
%        remaining optional arguments are passed through to fit().
%
% :Outputs:
%
%   **obj:**
%        The glm_map with the new basis set, a rebuilt design, stale fields
%        cleared, and either re-fit (if data given) or with refreshed design
%        diagnostics.
%
% :Examples:
% ::
%
%     d = fmri_glm_design_matrix(2, 'nscan', 200, 'units', 'secs', ...
%             'onsets', {[10 40 70]' [25 55 85]'}, 'condition_names', {'A' 'B'});
%     g = glm_map(d); g = build_design(g);
%     [xBF_hires, ~] = fmri_spline_basis(2, 'length', 20, 'nbasis', 4, 'order', 3);
%     g = replace_basis_set(g, 1, xBF_hires);                 % condition A -> spline
%     g = replace_basis_set(g, 1, xBF_hires, 'data', bold);   % and re-fit
%
% :See also:
%   - fmri_glm_design_matrix.replace_basis_set, build_design, fit, diagnostics
%
% ..
%    2026 - Initial implementation.
% ..

% -------------------------------------------------------------------------
% Validate
% -------------------------------------------------------------------------
if isempty(obj.design) || ~isa(obj.design, 'fmri_glm_design_matrix')
    error('glm_map:NoEventDesign', ...
        ['replace_basis_set requires an event/1st-level glm_map with an ' ...
         'fmri_glm_design_matrix in obj.design.']);
end

% Optional data for re-fitting (and pass-through fit options)
data = [];
wh = find(strcmpi(varargin, 'data'));
if ~isempty(wh)
    data = varargin{wh(1) + 1};
    varargin(wh(1):wh(1) + 1) = [];
end

had_contrasts = ~isempty(obj.contrasts);

% -------------------------------------------------------------------------
% Replace the basis set in the wrapped design and rebuild X
% -------------------------------------------------------------------------
obj.design = replace_basis_set(obj.design, condition_num, xBF_hires);
obj.design = build(obj.design);
obj.level  = 1;

% -------------------------------------------------------------------------
% Clear fields left over from the previous basis set that no longer apply.
% Changing the basis set changes the number and meaning of design columns,
% so fitted maps, contrasts, the fit's input structs, and diagnostics are
% all invalidated.
% -------------------------------------------------------------------------
obj.betas              = [];
obj.t                  = [];
obj.contrast_estimates = [];
obj.contrast_t         = [];
obj.df                 = [];
obj.sigma              = [];
obj.residuals          = [];
obj.dfe                = [];

obj.contrasts          = [];     % defined over the old regressors
obj.contrast_names     = {};
obj.contrast_summary_table = table();

obj.input_parameters     = struct();
obj.input_image_metadata = struct();
obj.fit_parameters       = struct();
obj.diagnostics          = struct();

% Restore the full (empty) schema for the nested structs
obj = validate_object(obj);

if had_contrasts
    obj.warnings{end + 1} = ['replace_basis_set: contrasts were cleared because the ' ...
        'design columns changed with the new basis set. Re-add contrasts with add_contrasts.'];
end

obj.history{end + 1} = sprintf('replace_basis_set: condition %d -> %s [%d x %d design]', ...
    condition_num, local_bf_name(xBF_hires), size(obj.X, 1), size(obj.X, 2));

% -------------------------------------------------------------------------
% Re-fit if data supplied; otherwise refresh design diagnostics
% -------------------------------------------------------------------------
if ~isempty(data)
    obj = fit(obj, data, varargin{:});
elseif ~isempty(obj.X)
    obj = run_diagnostics(obj, 'noverbose');
end

end % replace_basis_set


% =========================================================================
function nm = local_bf_name(xBF)
if isstruct(xBF) && isfield(xBF, 'name') && ~isempty(xBF.name)
    nm = xBF.name;
else
    nm = 'custom basis set';
end
end
