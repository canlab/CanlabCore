function obj = from_metadata_distance(meta, col, varargin)
% rsm.from_metadata_distance  Continuous-distance model RDM from a metadata column.
%
% Builds a model RDM where entry (i, j) = distance(v(i), v(j)) for a numeric
% metadata column v. Lets users model ordinal/continuous predictors (e.g.
% session distance, time since baseline, run distance) as continuous terms in
% the LME design rather than collapsing them to same-vs-different.
%
% Usage
% -----
%   M = rsm.from_metadata_distance(meta_table, 'session_number')
%   M = rsm.from_metadata_distance(meta_table, 'session_number', ...
%           'metric', 'abs_diff')
%   M = rsm.from_metadata_distance(numeric_vec, '', 'metric', 'squared_diff')
%
% Inputs
% ------
%   meta     Either a table (with `col` naming a numeric column) or a
%            numeric vector of length k.
%   col      Column name (when meta is a table). Ignored when meta is a vector.
%
%   varargin name-value pairs:
%       'metric'  one of {'abs_diff','squared_diff','log_diff','custom'}.
%                 Default: 'abs_diff'.
%       'fcn'     function handle (a, b) -> distance scalar. Required when
%                 'metric' is 'custom'. Default: [].
%       'labels'  {k x 1} cellstr — condition labels.
%       'name'    char — short name for this model RDM (stored in .source).
%
% Output
% ------
%   obj   rsm object with is_dissimilarity=true, metric='distance'

p = inputParser;
p.addParameter('metric', 'abs_diff', @(x) ischar(x) || isstring(x));
p.addParameter('fcn',    [],         @(x) isempty(x) || isa(x, 'function_handle'));
p.addParameter('labels', {},         @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
p.addParameter('name',   '',         @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

if istable(meta)
    if isempty(col) || ~ischar(col) && ~isstring(col)
        error('rsm:from_metadata_distance:badCol', ...
            'When meta is a table, you must pass a column name as the 2nd argument.');
    end
    v = meta.(char(col));
    if isempty(p.Results.labels)
        % Pre-assign so we don't overwrite with default later
        labels_default = cellstr(string(v));
    else
        labels_default = {};
    end
else
    v = meta;
    labels_default = {};
end

if ~isnumeric(v)
    try
        v = double(v);
    catch
        error('rsm:from_metadata_distance:nonNumeric', ...
            'metadata column must be numeric (or coercible to double).');
    end
end

v = v(:);
k = numel(v);

metric_name = lower(char(p.Results.metric));
switch metric_name
    case 'abs_diff'
        M = abs(v - v');
    case 'squared_diff'
        M = (v - v').^2;
    case 'log_diff'
        % log|a - b + eps|; pad tiny so identical values give 0 distance not -Inf
        D = abs(v - v');
        M = log1p(D);
    case 'custom'
        if isempty(p.Results.fcn)
            error('rsm:from_metadata_distance:noFcn', ...
                'metric=''custom'' requires a function handle via ''fcn''.');
        end
        M = zeros(k);
        for i = 1:k
            for j = i+1:k
                M(i, j) = p.Results.fcn(v(i), v(j));
                M(j, i) = M(i, j);
            end
        end
    otherwise
        error('rsm:from_metadata_distance:badMetric', ...
            'metric must be one of {abs_diff, squared_diff, log_diff, custom}; got %s.', metric_name);
end

M(1:k+1:end) = 0;

labels = cellstr(p.Results.labels);
if isempty(labels), labels = labels_default; end

name = char(p.Results.name);
if isempty(name)
    if istable(meta), name = char(col); else, name = 'distance'; end
end

obj = rsm(M, ...
    'is_dissimilarity', true, ...
    'metric',           'distance', ...
    'labels',           labels, ...
    'level',            'model_stack', ...
    'source',           ['design:' name '(' metric_name ')']);

end
