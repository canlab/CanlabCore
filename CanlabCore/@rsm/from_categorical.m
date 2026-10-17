function obj = from_categorical(v, varargin)
% rsm.from_categorical  Same-vs-different binary RDM from a categorical vector.
%
% Builds a model RDM where entry (i, j) = 0 iff v(i) == v(j), else 1. The
% diagonal is 0. This is the "SameX" predictor used as a building block for
% Phase 3 LME design matrices and for hypothesis RDMs (Bodysite, Condition,
% etc.) in the Sun et al. workflow.
%
% Wraps rsa.rdm.categoricalRDM when available, falls back to the equivalent
% `double(v(:) ~= v(:)')` otherwise.
%
% Usage
% -----
%   M = rsm.from_categorical(v)
%   M = rsm.from_categorical(v, 'labels', labels, 'name', 'Bodysite')
%
%   % Multi-column form: pass a table + list of column names, get an array
%   % of rsm objects, one per column.
%   M = rsm.from_categorical(meta_table, ...
%         {'subject_id','session_number','condition','bodysite'})
%
% Inputs
% ------
%   v        Either:
%              - vector (numeric / cell / categorical / string) of length k, OR
%              - table whose rows describe k conditions (used with column list)
%
%   varargin name-value pairs:
%       'columns'  cellstr — column names to pull from v (when v is a table).
%                  If omitted and v is a table with a 2nd positional argument
%                  that is cellstr, that arg is used as columns.
%       'labels'   {k x 1} cellstr — condition labels for the resulting RDM.
%       'name'     char    — short name for this model RDM (stored in .source).
%
% Output
% ------
%   obj   rsm object (or [1 x p] array when multiple columns requested)

% --------------------------------------------------------------
% Multi-column dispatch when first input is a table
% --------------------------------------------------------------
if istable(v)
    cols = [];
    if ~isempty(varargin) && (iscellstr(varargin{1}) || isstring(varargin{1})) %#ok<ISCLSTR>
        cols     = cellstr(varargin{1});
        varargin = varargin(2:end);
    elseif ~isempty(varargin) && ischar(varargin{1}) && ...
            ismember(varargin{1}, v.Properties.VariableNames)
        % Single char column name (e.g. from_categorical(meta, 'condition'))
        cols     = varargin(1);
        varargin = varargin(2:end);
    end
    p = inputParser;
    p.addParameter('columns', cols, @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
    p.addParameter('labels', {},   @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
    p.addParameter('name', '',     @(x) ischar(x) || isstring(x));
    p.parse(varargin{:});
    cols = cellstr(p.Results.columns);

    if isempty(cols)
        error('rsm:from_categorical:noColumns', ...
            'When passing a table, you must also specify which column(s) to use.');
    end

    obj = rsm.empty;
    for c = 1:numel(cols)
        sub_name = cols{c};
        obj(end+1) = rsm.from_categorical(v.(sub_name), ...
            'labels', p.Results.labels, ...
            'name',   sub_name); %#ok<AGROW>
    end
    return
end

% --------------------------------------------------------------
% Single-vector path
% --------------------------------------------------------------
p = inputParser;
p.addParameter('labels', {}, @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
p.addParameter('name', '',   @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

v = v(:);
k = numel(v);

caps = probe_rsatoolbox();

% Build the binary same-vs-different matrix
if caps.rdm_categoricalRDM
    try
        if ~isempty(which('rsa.rdm.categoricalRDM'))
            M = rsa.rdm.categoricalRDM(v);
        else
            M = categoricalRDM(v);
        end
        % rsa.rdm.categoricalRDM may return packed structures depending on
        % version; coerce to a plain k x k binary matrix.
        if isstruct(M) && isfield(M, 'RDM'), M = M.RDM; end
        M = double(M);
        % Some rsatoolbox versions return a "categorical crossing count" rather
        % than a binary same-diff. Re-binarize defensively.
        M = double(M ~= 0);
    catch
        M = local_categorical_rdm(v);
    end
else
    M = local_categorical_rdm(v);
end

% Diagonal is always 0 (same-as-self)
M(1:k+1:end) = 0;

labels = cellstr(p.Results.labels);
if isempty(labels)
    % Default labels: stringified values of v
    if iscell(v)
        labels = cellfun(@(x) char(string(x)), v, 'UniformOutput', false);
    else
        labels = cellstr(string(v));
    end
end

name = char(p.Results.name);
if isempty(name), name = 'categorical'; end

obj = rsm(M, ...
    'is_dissimilarity', true, ...
    'metric',           'categorical', ...
    'labels',           labels, ...
    'level',            'model_stack', ...
    'source',           ['design:' name]);

end


% =========================================================================
function M = local_categorical_rdm(v)
% Stock fallback: 1 where v(i) ~= v(j), else 0; diagonal = 0.
    if iscategorical(v)
        % Convert to double codes for fast comparison
        u = double(v);
    elseif iscell(v)
        % Map unique values to integers
        [~, ~, u] = unique(v, 'stable');
    elseif isstring(v)
        [~, ~, u] = unique(v, 'stable');
    else
        u = double(v);
    end
    u = u(:);
    M = double(u ~= u');
end
