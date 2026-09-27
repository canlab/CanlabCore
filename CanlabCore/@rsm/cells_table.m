function T = cells_table(obj, spec, varargin)
% cells_table  Pull multiple grouping pairs into a single per-replicate table.
%
% Wraps R.cells() and stacks the results into a table with one column per
% requested grouping (or grouping pair). Rows correspond to replicates
% (typically subjects), so the table is directly consumable by ttest,
% barplot_columns, fitlme, etc.
%
% Usage
% -----
%   % Shorthand: a list of grouping names — each is treated as within-group
%   %   ('hot' -> R.cells('hot','hot'), etc.)
%   T = R.cells_table({'hot','warm','imagine'})
%
%   % Mixed: scalar (within) and 1x2 cells (between)
%   T = R.cells_table({ ...
%       'hot', ...                                   % within hot
%       'warm', ...                                  % within warm
%       {'hot','warm'}, ...                          % between hot and warm
%       {'hot','imagine'}, ...                       % between hot and imagine
%       struct('name','HW', 'a','hot', 'b','warm') ...  % named entry
%   })
%
%   % With explicit names
%   T = R.cells_table({'hot','warm','imagine'}, 'names',{'Hot','Warm','Imagine'})
%
% Inputs
% ------
%   spec     cell array; each entry is one of:
%              char/string         -- within-group (uses the name on both sides)
%              {a, b}              -- between-group (1x2 cell of names/indices)
%              struct(.a,.b,.name) -- explicit with name
%   varargin name-value (forwarded to R.cells):
%              'transform'   'auto' | 'fisherz' | 'none'   (default 'auto')
%              'reduction'   'mean' | 'median' | 'sum'      (default 'mean')
%              'names'       cellstr of column names overriding auto-derived ones
%
% Output
% ------
%   T  table with one column per spec entry; height = N_replicates.

if builtin('numel', obj) > 1
    error('rsm:cells_table:nonScalar', 'cells_table() expects a scalar rsm.');
end

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('names', {}, @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
p.parse(varargin{:});
override_names = cellstr(p.Results.names);

% Forward only the non-'names' params to cells()
forward = remove_param(varargin, 'names');

T = table;
n_entries = numel(spec);
auto_names = cell(n_entries, 1);

for i = 1:n_entries
    entry = spec{i};

    if ischar(entry) || isstring(entry)
        a = char(entry); b = a;
        nm = char(entry);
    elseif iscell(entry) && numel(entry) == 2
        a = entry{1}; b = entry{2};
        nm = sprintf('%s_%s', stringify(a), stringify(b));
    elseif isstruct(entry)
        a = entry.a; b = entry.b;
        if isfield(entry, 'name') && ~isempty(entry.name), nm = char(entry.name);
        else, nm = sprintf('%s_%s', stringify(a), stringify(b));
        end
    else
        error('rsm:cells_table:badSpec', ...
            'spec entry %d must be a name (within), {a, b} cell, or struct.', i);
    end

    v = obj.cells(a, b, forward{:});

    auto_names{i} = sanitize_varname(nm);
    T.(auto_names{i}) = v;
end

% Apply user-supplied names if provided
if ~isempty(override_names)
    if numel(override_names) ~= n_entries
        error('rsm:cells_table:badNames', ...
            'numel(names) = %d but spec has %d entries.', numel(override_names), n_entries);
    end
    T.Properties.VariableNames = cellfun(@sanitize_varname, override_names, 'UniformOutput', false);
end

end


function s = stringify(x)
if ischar(x) || isstring(x), s = char(x);
elseif isnumeric(x), s = ['idx' strrep(num2str(x(:)'), ' ','')];
elseif islogical(x), s = 'mask';
else, s = 'group';
end
end


function v = sanitize_varname(name)
% Replace any non-word chars with underscore so it's a valid table column name
v = regexprep(char(name), '[^A-Za-z0-9_]', '_');
if isempty(v) || ~isletter(v(1)), v = ['x' v]; end
end


function out = remove_param(args, name)
% Strip a name-value pair from a varargin cell
out = args;
idx = find(cellfun(@(x) (ischar(x) || isstring(x)) && strcmpi(char(x), name), out));
if ~isempty(idx) && idx < numel(out)
    out([idx, idx+1]) = [];
end
end
