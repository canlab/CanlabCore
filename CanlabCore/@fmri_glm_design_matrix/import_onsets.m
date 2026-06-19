function obj = import_onsets(obj, source, varargin)
% Import event onsets/durations/condition names (and parametric modulators)
% into an fmri_glm_design_matrix object. Two input styles are supported.
%
% (1) Tabular / FSL-style -- a file (.csv/.txt/.xlsx), or a MATLAB table, with
%     one row per event and columns for onset, (optional) duration, and a
%     condition label. The label column may be a string name OR an integer
%     event-type code (condition and event type mean the same thing). Events
%     are grouped by label into conditions (one Sess.U entry each).
%
% (2) SPM-style -- cell arrays with one cell per condition: a cell of onset
%     vectors, and (optionally) matching cells of durations and parametric
%     modulator values.
%
% :Usage:
% ::
%
%     % Tabular / FSL-style
%     obj = import_onsets(obj, 'events.csv')
%     obj = import_onsets(obj, events_table, 'onset_col','Onset', 'name_col','EV')
%
%     % SPM-style cell arrays (one cell per condition / event type)
%     obj = import_onsets(obj, onsets_cell)
%     obj = import_onsets(obj, onsets_cell, durations_cell)
%     obj = import_onsets(obj, onsets_cell, durations_cell, pmods_cell, ...
%                         'names', {'A','B'}, 'pm_names', {'rt',''})
%
% :Inputs:
%
%   **obj:**  an fmri_glm_design_matrix object.
%
%   **source:**
%        A filename, a MATLAB table (tabular style), OR a cell array of onset
%        vectors, one cell per condition (SPM style).
%
% :Optional Inputs (tabular style):
%
%   **'onset_col' / 'dur_col' / 'name_col', columnname:**
%        Override the column used for onsets / durations / condition labels.
%        Defaults are matched case-insensitively against common names
%        (onset/duration/name/condition/trial_type/event_type/value/code).
%
% :Optional Inputs (SPM style):
%
%   **second/third positional cells:**
%        durations (cell, one per condition) and parametric-modulator values
%        (cell, one per condition); pass [] to skip either.
%
%   **'names', {...}:**      condition names, one per condition.
%   **'pm_names', {...}:**   parametric-modulator names, one per condition.
%
% :Outputs:
%
%   **obj:**  with obj.Sess(1).U populated (one entry per condition).
%
% :Examples:
% ::
%
%     d = fmri_glm_design_matrix(2, 'nscan', 200, 'units', 'secs');
%     d = import_onsets(d, 'events.csv');                 % FSL/tabular
%     d = import_onsets(d, {[10 40]' [25 55]'}, {4 4});   % SPM-style
%     d = build(d);
%
% :See also: fmri_glm_design_matrix.add, Add_Event_Info, readtable
%
% ..
%    2026 - Reimplemented (previous version was a non-functional stub) and
%    extended to FSL/tabular files with event-type codes and SPM-style cells.
% ..

if iscell(source)
    obj = local_import_spm_cells(obj, source, varargin{:});
elseif istable(source) || ischar(source) || isstring(source)
    obj = local_import_table(obj, source, varargin{:});
else
    error('fmri_glm_design_matrix:import_onsets', ...
        'source must be a filename, a table, or a cell array of onset vectors.');
end

end % import_onsets


% =========================================================================
% SPM-style cell arrays
% =========================================================================
function obj = local_import_spm_cells(obj, onsets, varargin)

% Leading positional cell/[] args are durations then pmods
durations = {}; pmods = {};
k = 0;
while ~isempty(varargin) && (iscell(varargin{1}) || isempty(varargin{1})) && ~ischar(varargin{1})
    k = k + 1;
    if k == 1, durations = varargin{1}; elseif k == 2, pmods = varargin{1}; else, break, end
    varargin(1) = [];
end

p = inputParser;
p.addParameter('names', {}, @iscell);
p.addParameter('pm_names', {}, @iscell);
p.parse(varargin{:});
names = p.Results.names; pm_names = p.Results.pm_names;

nconds = numel(onsets);

for i = 1:nconds

    obj.Sess(1).U(i).ons = onsets{i}(:);

    % name
    if i <= numel(names) && ~isempty(names{i}), obj.Sess(1).U(i).name = names{i};
    else, obj.Sess(1).U(i).name = sprintf('Cond%d', i);
    end

    % duration
    if i <= numel(durations) && ~isempty(durations{i})
        obj.Sess(1).U(i).dur = durations{i}(:);
    else
        obj.Sess(1).U(i).dur = zeros(numel(onsets{i}), 1);
    end

    % parametric modulator (single per condition)
    if i <= numel(pmods) && ~isempty(pmods{i})
        P = struct('name', '', 'P', pmods{i}(:), 'h', 1, 'dur', [], 'i', []);
        if i <= numel(pm_names) && ~isempty(pm_names{i}), P.name = pm_names{i};
        else, P.name = sprintf('pm%d', i);
        end
        obj.Sess(1).U(i).P = P;
    end
end

obj.history{end + 1} = sprintf('import_onsets: %d conditions from SPM-style cell arrays', nconds);

end


% =========================================================================
% Tabular / FSL-style
% =========================================================================
function obj = local_import_table(obj, source, varargin)

p = inputParser;
p.addParameter('onset_col', '', @(x) ischar(x) || isstring(x));
p.addParameter('dur_col',   '', @(x) ischar(x) || isstring(x));
p.addParameter('name_col',  '', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opt = p.Results;

if istable(source)
    T = source;
else
    if exist(char(source), 'file') ~= 2
        error('fmri_glm_design_matrix:import_onsets', 'File not found: %s', char(source));
    end
    T = readtable(char(source));
end

vn = T.Properties.VariableNames;

onset_col = local_pick_col(vn, opt.onset_col, {'onset', 'onsets', 'ons', 'time'});
name_col  = local_pick_col(vn, opt.name_col,  {'name', 'condition', 'cond', 'trial_type', ...
                                               'event_type', 'eventtype', 'type', 'value', 'code', 'ev'});
dur_col   = local_pick_col(vn, opt.dur_col,   {'duration', 'durations', 'dur'});

if isempty(onset_col)
    error('fmri_glm_design_matrix:import_onsets', ...
        'Could not find an onset column. Columns: %s. Use the ''onset_col'' option.', strjoin(vn, ', '));
end
if isempty(name_col)
    error('fmri_glm_design_matrix:import_onsets', ...
        'Could not find a condition/event-type column. Columns: %s. Use the ''name_col'' option.', strjoin(vn, ', '));
end

onsets = T.(onset_col);
labels = T.(name_col);

% Labels may be string names or integer event-type codes; normalize to a
% cellstr (codes become 'Cond<code>').
if isnumeric(labels)
    labels = arrayfun(@(v) sprintf('Cond%g', v), labels, 'UniformOutput', false);
elseif iscategorical(labels) || isstring(labels)
    labels = cellstr(labels);
elseif ~iscell(labels)
    labels = cellstr(string(labels));
end

if ~isempty(dur_col), durs = T.(dur_col); else, durs = zeros(size(onsets)); end

% Group by condition label (first-appearance order)
[uconds, ~, grp] = unique(labels, 'stable');

for c = 1:numel(uconds)
    wh = (grp == c);
    obj.Sess(1).U(c).name = uconds{c};
    obj.Sess(1).U(c).ons  = onsets(wh);
    obj.Sess(1).U(c).dur  = durs(wh);
end

obj.history{end + 1} = sprintf('import_onsets: %d conditions, %d events from %s', ...
    numel(uconds), numel(onsets), local_src_name(source));

end


% =========================================================================
function col = local_pick_col(varnames, override, candidates)
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


function s = local_src_name(source)
if istable(source), s = '(table)'; else, s = char(source); end
end
