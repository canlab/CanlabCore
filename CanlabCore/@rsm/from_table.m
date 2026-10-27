function obj = from_table(meta, spec, varargin)
% rsm.from_table  Build multiple model RDMs from a metadata table + spec.
%
% Convenience wrapper that dispatches each spec entry to from_categorical or
% from_metadata_distance based on a 'kind' field, returning an array of rsm
% objects.
%
% Usage
% -----
%   spec = {
%       struct('col','condition',      'kind','categorical');
%       struct('col','bodysite',       'kind','categorical');
%       struct('col','session_number', 'kind','distance', 'metric','abs_diff');
%   };
%   M = rsm.from_table(meta_table, spec);
%
% Shorthand (auto-detects: numeric column -> distance; everything else ->
% categorical):
%   M = rsm.from_table(meta_table, {'condition','bodysite','session_number'});
%
% Inputs
% ------
%   meta       table — one row per condition.
%   spec       Either:
%              - cell array of struct entries with fields:
%                  .col     column name in meta
%                  .kind    'categorical' or 'distance'
%                  .metric  (optional) metric for 'distance' kind
%                  .name    (optional) override the RDM's short name
%              - cellstr/string of column names (auto-detect kind)
%
%   varargin name-value pairs:
%       'labels'  {k x 1} cellstr — condition labels applied to all RDMs.
%
% Output
% ------
%   obj   [1 x p] array of rsm objects

p = inputParser;
p.addParameter('labels', {}, @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
p.parse(varargin{:});

% --------------------------------------------------------------
% Normalize spec into an array of structs with .col, .kind, .metric, .name
% --------------------------------------------------------------
if iscellstr(spec) || isstring(spec) %#ok<ISCLSTR>
    cols  = cellstr(spec);
    entries = cell(numel(cols), 1);
    for i = 1:numel(cols)
        c = cols{i};
        v = meta.(c);
        kind = 'categorical';
        if isnumeric(v) && ~all(v == round(v)) || ...
           (isnumeric(v) && numel(unique(v(:))) > 8 && all(v == round(v)))
            % Heuristic: numeric with non-integer values OR integer with many
            % levels -> treat as distance. Users wanting categorical instead
            % should use explicit struct spec.
            kind = 'distance';
        end
        entries{i} = struct('col', c, 'kind', kind, 'metric', 'abs_diff', 'name', c);
    end
elseif iscell(spec)
    entries = spec;
elseif isstruct(spec)
    entries = num2cell(spec);
else
    error('rsm:from_table:badSpec', ...
        'spec must be cellstr, cell-of-structs, or struct array.');
end

% --------------------------------------------------------------
% Dispatch
% --------------------------------------------------------------
obj = rsm.empty;
for i = 1:numel(entries)
    e = entries{i};
    if ~isfield(e, 'kind'), e.kind = 'categorical'; end
    if ~isfield(e, 'name') || isempty(e.name), e.name = e.col; end

    switch lower(e.kind)
        case 'categorical'
            obj(end+1) = rsm.from_categorical(meta.(e.col), ...
                'labels', p.Results.labels, ...
                'name',   e.name); %#ok<AGROW>
        case 'distance'
            metric_name = 'abs_diff';
            if isfield(e, 'metric') && ~isempty(e.metric), metric_name = e.metric; end
            obj(end+1) = rsm.from_metadata_distance(meta, e.col, ...
                'metric', metric_name, ...
                'labels', p.Results.labels, ...
                'name',   e.name); %#ok<AGROW>
        otherwise
            error('rsm:from_table:badKind', ...
                'spec.kind must be ''categorical'' or ''distance''; got %s.', e.kind);
    end
end

end
