function obj = validate_object(obj)
% Ensure a glm_map's nested structs expose their full set of valid fields.
%
% glm_map groups related outputs into nested structs (.input_parameters,
% .input_image_metadata, .diagnostics, and .fit_parameters). Depending on how
% an object was built (empty constructor, re-cast from an fmri_data.regress
% results struct, or partially populated), some of those structs may be empty
% or missing fields. validate_object fills in any missing fields with empty
% defaults so that every struct always exposes the complete schema, and
% normalizes the nested diagnostics.collinearity_report the same way.
%
% Existing field values are never overwritten; only absent fields are added.
% A struct property that is not a struct (e.g. left as []) is replaced by the
% empty-field template. The field templates are the single source of truth for
% the valid field names produced by fmri_data.regress and glm_map.diagnostics.
%
% :Usage:
% ::
%
%     obj = validate_object(obj)
%
% :Inputs:
%
%   **obj:**
%        A glm_map object.
%
% :Outputs:
%
%   **obj:**
%        The glm_map object with .input_parameters, .input_image_metadata,
%        .diagnostics (and .diagnostics.collinearity_report), and
%        .fit_parameters guaranteed to contain every valid field.
%
% :Examples:
% ::
%
%     g = validate_object(glm_map);          % empty object, full schema
%     fieldnames(g.diagnostics)              % all diagnostic fields present
%
% :See also:
%   - glm_map, glm_map.check_properties, fmri_data.regress
%
% ..
%    Programmers' notes:
%    2026 - Initial implementation. Templates below define the valid field
%    set for each nested struct; keep them in sync with the fields produced by
%    fmri_data.regress (out.input_parameters / out.input_image_metadata /
%    out.diagnostics) and glm_map.diagnostics.
% ..

obj.input_parameters     = local_fill(obj.input_parameters,     local_tmpl_input_parameters());
obj.input_image_metadata = local_fill(obj.input_image_metadata, local_tmpl_input_image_metadata());
obj.diagnostics          = local_fill(obj.diagnostics,          local_tmpl_diagnostics());
obj.fit_parameters       = local_fill(obj.fit_parameters,       local_tmpl_fit_parameters());

% Nested struct inside .diagnostics
obj.diagnostics.collinearity_report = local_fill(obj.diagnostics.collinearity_report, ...
    local_tmpl_collinearity_report());

end % validate_object


% =====================================================================
% Field-filling helper
% =====================================================================
function s = local_fill(s, tmpl)
% Add every field of tmpl that is missing from s (default = template value).
% Existing fields are preserved. A non-struct (or non-scalar struct) s is
% replaced by the template.
if ~isstruct(s) || ~isscalar(s)
    s = tmpl;
    return
end

fn = fieldnames(tmpl);
for i = 1:numel(fn)
    if ~isfield(s, fn{i})
        s.(fn{i}) = tmpl.(fn{i});
    end
end
end


% =====================================================================
% Field templates (valid field sets). Defaults are [] = "not set / not
% computed". Keep in sync with fmri_data.regress and glm_map.diagnostics.
% =====================================================================
function t = local_tmpl_input_parameters()
t = struct( ...
    'brain_is_predictor',            [], ...
    'do_robust',                     [], ...
    'grandmeanscale',                [], ...
    'do_intercept',                  [], ...
    'do_resid',                      [], ...
    'doverbose',                     [], ...
    'do_display',                    [], ...
    'covdat',                        [], ...
    'initial_statistical_threshold', []);
end


function t = local_tmpl_input_image_metadata()
t = struct( ...
    'source_notes', [], ...
    'history',      [], ...
    'image_names',  [], ...
    'fullpath',     []);
end


function t = local_tmpl_diagnostics()
t = struct( ...
    'Variance_inflation_factors',          [], ...
    'Contrast_variance_inflation_factors', [], ...
    'Leverages',                           [], ...
    'Cooks_distance',                      [], ...
    'condition_number',                    [], ...
    'rank_deficient',                      [], ...
    'collinearity_report',                 [], ...
    'vif_threshold',                       []);
end


function t = local_tmpl_collinearity_report()
t = struct( ...
    'vif_threshold',          [], ...
    'high_vif_columns',       [], ...
    'duplicate_column_pairs', [], ...
    'high_correlation_pairs', []);
end


function t = local_tmpl_fit_parameters()
t = struct( ...
    'robust',        [], ...
    'ar_order',      [], ...
    'is_timeseries', [], ...
    'pthresh',       [], ...
    'thresh_type',   [], ...
    'do_resid',      []);
end
