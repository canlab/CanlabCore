function obj = import_onsets(obj, source, varargin)
% Import event onsets into a glm_map's wrapped fmri_glm_design_matrix design.
%
% Passes its inputs through to fmri_glm_design_matrix.import_onsets (which
% accepts FSL/tabular files or tables, and SPM-style cell arrays), then builds
% the design matrix if enough information is available. Event regressors are
% flagged as of interest automatically (they form the design's H partition);
% any covariates of no interest become nuisance regressors. If the glm_map has
% no design yet, you can bootstrap one by passing 'TR' and 'nscan'; otherwise
% the method explains what is needed and how to attach it.
%
% :Usage:
% ::
%
%     g = import_onsets(g, source, ...)                       % design already attached
%     g = import_onsets(g, source, 'TR', 2, 'nscan', 200)     % bootstrap the design
%
% :Inputs:
%
%   **obj:**     a glm_map object.
%   **source:**  a filename / table (FSL/tabular) or a cell of onset vectors
%                (SPM-style); see fmri_glm_design_matrix.import_onsets.
%
% :Optional Inputs:
%
%   **'TR', value / 'nscan', value / 'units', 'secs'|'scans':**
%        Used to create (or update) the wrapped design when bootstrapping.
%        Remaining options are passed through to the design's import_onsets.
%
% :Outputs:
%
%   **obj:**  glm_map with the design imported and (if possible) built; level
%             set to 1. obj.wh_interest then flags the event regressors.
%
% :Examples:
% ::
%
%     g = import_onsets(glm_map, 'events.csv', 'TR', 2, 'nscan', 200);
%     g = run_diagnostics(g);
%
% :See also:
%   - fmri_glm_design_matrix.import_onsets, build_design, import_SPM
%
% ..
%    2026 - Initial implementation.
% ..

% -------------------------------------------------------------------------
% Pull out design-bootstrap options (TR / nscan / units); pass the rest on
% -------------------------------------------------------------------------
[TR, nscan, units, rest] = local_extract_design_opts(varargin);

% -------------------------------------------------------------------------
% Ensure a design object exists
% -------------------------------------------------------------------------
if isempty(obj.design) || ~isa(obj.design, 'fmri_glm_design_matrix')
    if ~isempty(TR) && ~isempty(nscan)
        d = fmri_glm_design_matrix(TR, 'nscan', nscan);
        if ~isempty(units), d = add(d, 'units', units); end
        obj.design = d;
    else
        error('glm_map:NoDesignForImport', '%s', local_design_guidance());
    end
else
    if ~isempty(TR),    obj.design.TR = TR; obj.design.xY.RT = TR; end
    if ~isempty(nscan), obj.design.nscan = nscan; end
    if ~isempty(units), obj.design = add(obj.design, 'units', units); end
end

% -------------------------------------------------------------------------
% Delegate the import to the wrapped design object
% -------------------------------------------------------------------------
obj.design = import_onsets(obj.design, source, rest{:});
obj.level  = 1;

% -------------------------------------------------------------------------
% Build if we have enough information; otherwise direct the user
% -------------------------------------------------------------------------
if local_can_build(obj.design)
    obj.design = build(obj.design);
    obj = validate_object(obj);
    obj.history{end + 1} = sprintf('import_onsets: imported and built design [%d x %d]', ...
        size(obj.X, 1), size(obj.X, 2));
else
    obj.history{end + 1} = 'import_onsets: imported onsets; design not built (missing info)';
    fprintf('%s\n', local_missing_info(obj.design));
end

end % import_onsets


% =========================================================================
% Local helpers
% =========================================================================
function [TR, nscan, units, rest] = local_extract_design_opts(args)
[TR, nscan, units] = deal([]);
rest = args;
for key = {'TR', 'nscan', 'units'}
    wh = find(strcmpi(rest, key{1}), 1);
    if ~isempty(wh) && wh < numel(rest)
        switch lower(key{1})
            case 'tr',    TR = rest{wh + 1};
            case 'nscan', nscan = rest{wh + 1};
            case 'units', units = rest{wh + 1};
        end
        rest(wh:wh + 1) = [];
    end
end
end


function tf = local_can_build(d)
% Buildable if TR is set, nscan is set, and at least one condition has onsets.
tf = false;
if isempty(d) || ~isa(d, 'fmri_glm_design_matrix'), return, end
if isempty(d.nscan) || any(isnan(d.TR)), return, end
if isempty(d.Sess) || ~isfield(d.Sess(1), 'U') || isempty(d.Sess(1).U), return, end
if ~isfield(d.Sess(1).U(1), 'ons') || isempty(d.Sess(1).U(1).ons), return, end
tf = true;
end


function s = local_missing_info(d)
miss = {};
if isempty(d.nscan) || any(isnan(d.TR))
    if any(isnan(d.TR)),   miss{end + 1} = 'TR (repetition time)'; end
    if isempty(d.nscan),   miss{end + 1} = 'nscan (scans per session)'; end
end
if isempty(d.Sess) || ~isfield(d.Sess(1), 'U') || isempty(d.Sess(1).U) ...
        || ~isfield(d.Sess(1).U(1), 'ons') || isempty(d.Sess(1).U(1).ons)
    miss{end + 1} = 'event onsets';
end
s = sprintf([ ...
    '  import_onsets: design imported but not built; still missing: %s.\n' ...
    '  Attach the missing info and build, e.g.:\n' ...
    '    g.design.nscan = nscan;   %% number of scans per session\n' ...
    '    g.design.TR    = TR;      %% repetition time (s)\n' ...
    '    g = build_design(g);'], strjoin(miss, ', '));
end


function s = local_design_guidance()
s = sprintf([ ...
    'glm_map has no design to import onsets into.\n' ...
    'Either bootstrap one by passing TR and nscan:\n' ...
    '    g = import_onsets(g, source, ''TR'', 2, ''nscan'', 200, ''units'', ''secs'');\n' ...
    'or attach an fmri_glm_design_matrix first:\n' ...
    '    g.design = fmri_glm_design_matrix(2, ''nscan'', 200, ''units'', ''secs'');\n' ...
    '    g = import_onsets(g, source);']);
end
