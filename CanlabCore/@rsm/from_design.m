function obj = from_design(X, varargin)
% rsm.from_design  Model RDM(s) from a numeric/binary design matrix.
%
% Reproduces the design-to-RDM logic that @fmri_data/rsa_regression.m embeds
% inline (its lines 95–100):
%       modelRDM(:, i) = pdist(design(:, i), 'seuclidean')
%       modelRDM(:, i) = 100000 * modelRDM(:, i) / sum(modelRDM(:, i))
%
% Each column of X is treated as one regressor and produces one model RDM.
%
% Usage
% -----
%   M = rsm.from_design(X)
%   M = rsm.from_design(X, 'names', {'hot','warm','imagine'}, ...
%                          'metric', 'seuclidean', ...
%                          'normalize', true)
%
% Inputs
% ------
%   X        [k x p] numeric design matrix. k = number of conditions/observations,
%            p = number of regressors. Each column becomes one model RDM.
%
%   varargin name-value pairs:
%       'names'      {1 x p} cellstr — names for each model RDM (default: regressor_i).
%       'metric'     pdist metric, default 'seuclidean'. Common alternatives:
%                    'euclidean', 'cityblock', 'hamming'.
%       'normalize'  logical, default true. Match rsa_regression's
%                    100000 * d / sum(d) rescaling.
%       'labels'     {k x 1} cellstr — condition (row/col) labels.
%
% Output
% ------
%   obj   [1 x p] array of rsm objects (or scalar if p == 1)

p = inputParser;
p.addParameter('names',     {},          @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
p.addParameter('metric',    'seuclidean',@(x) ischar(x) || isstring(x));
p.addParameter('normalize', true,        @(x) islogical(x) || isnumeric(x));
p.addParameter('labels',    {},          @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
p.parse(varargin{:});

if ~isnumeric(X) && ~islogical(X)
    error('rsm:from_design:nonNumeric', 'X must be a numeric or logical [k x p] design matrix.');
end
X = double(X);

[k, q] = size(X);

names = cellstr(p.Results.names);
if isempty(names)
    names = arrayfun(@(i) sprintf('regressor_%d', i), 1:q, 'UniformOutput', false);
elseif numel(names) ~= q
    error('rsm:from_design:badNames', ...
        'numel(names) = %d but X has %d columns.', numel(names), q);
end

labels = cellstr(p.Results.labels);
if isempty(labels)
    labels = arrayfun(@(i) sprintf('cond_%d', i), 1:k, 'UniformOutput', false)';
end

metric_name = char(p.Results.metric);
normalize   = logical(p.Results.normalize);

obj = rsm.empty;
for i = 1:q
    d = pdist(X(:, i), metric_name);   % 1 x [k*(k-1)/2]
    if normalize && sum(d) > 0
        d = 100000 * d / sum(d);
    end
    M = squareform(d);
    M(1:k+1:end) = 0;

    obj(end+1) = rsm(M, ...                                              %#ok<AGROW>
        'is_dissimilarity', true, ...
        'metric',           'design', ...
        'labels',           labels, ...
        'level',            'model_stack', ...
        'source',           ['design:' names{i} '(' metric_name ')']);
end

if q == 1, obj = obj(1); end

end
