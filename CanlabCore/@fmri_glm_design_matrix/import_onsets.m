function obj = import_onsets(obj, design_file, varargin)
% Import event onsets/durations/condition names from a table file and attach
% them to an fmri_glm_design_matrix object (single session).
%
% :Usage:
% ::
%
%     obj = import_onsets(obj, design_file)
%     obj = import_onsets(obj, design_file, 'onset_col','Onset', 'dur_col','Dur', 'name_col','Trial')
%
% :Inputs:
%
%   **obj:**
%        An fmri_glm_design_matrix object.
%
%   **design_file:**
%        Path to a .csv / .txt / .xlsx file, or a MATLAB table. It must have
%        one row per event and (by default) columns named onset, duration,
%        and name (case-insensitive; "condition" or "trial_type" are also
%        accepted for the name column). Events are grouped by name into
%        conditions; onsets and durations are read in the file's units (the
%        same units, secs or scans, that the object's basis set expects).
%
% :Optional Inputs:
%
%   **'onset_col' / 'dur_col' / 'name_col', columnname:**
%        Override the column name used for onsets / durations / condition
%        labels. dur_col may be omitted from the file (durations default 0).
%
% :Outputs:
%
%   **obj:**
%        The object with obj.Sess(1).U populated (one entry per condition).
%
% :Examples:
% ::
%
%     d = fmri_glm_design_matrix(2, 'nscan', 200, 'units', 'secs');
%     d = import_onsets(d, 'events.csv');
%     d = build(d);
%
% :See also: fmri_glm_design_matrix.add, Add_Event_Info, readtable
%
% ..
%    2026 - Reimplemented (previous version was a non-functional stub).
% ..

% -------------------------------------------------------------------------
% Parse options
% -------------------------------------------------------------------------
p = inputParser;
p.addParameter('onset_col', '', @(x) ischar(x) || isstring(x));
p.addParameter('dur_col',   '', @(x) ischar(x) || isstring(x));
p.addParameter('name_col',  '', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opt = p.Results;

% -------------------------------------------------------------------------
% Read the table
% -------------------------------------------------------------------------
if istable(design_file)
    T = design_file;
else
    if ~ischar(design_file) && ~isstring(design_file)
        error('fmri_glm_design_matrix:import_onsets', 'design_file must be a filename or a table.');
    end
    if exist(design_file, 'file') ~= 2
        error('fmri_glm_design_matrix:import_onsets', 'File not found: %s', char(design_file));
    end
    T = readtable(design_file);
end

vn = T.Properties.VariableNames;

onset_col = local_pick_col(vn, opt.onset_col, {'onset', 'onsets', 'ons', 'time'});
name_col  = local_pick_col(vn, opt.name_col,  {'name', 'condition', 'cond', 'trial_type', 'type'});
dur_col   = local_pick_col(vn, opt.dur_col,   {'duration', 'durations', 'dur'});  % may be ''

if isempty(onset_col)
    error('fmri_glm_design_matrix:import_onsets', ...
        'Could not find an onset column. Columns: %s. Use the ''onset_col'' option.', strjoin(vn, ', '));
end
if isempty(name_col)
    error('fmri_glm_design_matrix:import_onsets', ...
        'Could not find a condition-name column. Columns: %s. Use the ''name_col'' option.', strjoin(vn, ', '));
end

onsets = T.(onset_col);
labels = T.(name_col);
if ~iscell(labels), labels = cellstr(string(labels)); end

if ~isempty(dur_col), durs = T.(dur_col); else, durs = zeros(size(onsets)); end

% -------------------------------------------------------------------------
% Group by condition name (preserving first-appearance order)
% -------------------------------------------------------------------------
[uconds, ~, grp] = unique(labels, 'stable');

for c = 1:numel(uconds)
    wh = (grp == c);
    obj.Sess(1).U(c).name = uconds{c};
    obj.Sess(1).U(c).ons  = onsets(wh);
    obj.Sess(1).U(c).dur  = durs(wh);
end

obj.history{end + 1} = sprintf('import_onsets: %d conditions, %d events from %s', ...
    numel(uconds), numel(onsets), local_src_name(design_file));

end % import_onsets


% =========================================================================
function col = local_pick_col(varnames, override, candidates)
% Resolve a column name: use override if given, else first case-insensitive
% match among candidates, else ''.
col = '';
if ~isempty(override)
    wh = find(strcmpi(varnames, override), 1);
    if isempty(wh)
        error('fmri_glm_design_matrix:import_onsets', 'Column ''%s'' not found.', char(override));
    end
    col = varnames{wh};
    return
end
for k = 1:numel(candidates)
    wh = find(strcmpi(varnames, candidates{k}), 1);
    if ~isempty(wh), col = varnames{wh}; return, end
end
end


function s = local_src_name(design_file)
if istable(design_file), s = '(table)'; else, s = char(design_file); end
end
