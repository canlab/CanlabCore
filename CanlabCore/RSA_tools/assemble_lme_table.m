function [tbl, info] = assemble_lme_table(dat, varargin)
% assemble_lme_table  Build a long-format table for LME modeling of RSA data.
%
% Implements the contract from `RSA_tools/RSA_Phase3_LME_Design.md` §3:
% one row per (i, j) upper-triangle pair from a single subject's image-level
% RSM. Columns include Y (similarity, Fisher-z by default), Subject
% (categorical grouping for random effects), and binary "Same<predictor>"
% columns for each metadata predictor + element-wise AND interaction columns.
%
% This is the engine behind @fmri_data/rsa_lme and @fmri_data/rsa_lm. It is
% also exported as a free function so users can build tables and run their
% own fitlme / fitlm calls.
%
% Usage
% -----
%   [tbl, info] = assemble_lme_table(dat, ...
%       'predictors',   {'condition','bodysite','session_number'}, ...
%       'interactions', {{'condition','bodysite'}}, ...
%       'subject_var',  'sub')
%
% Inputs
% ------
%   dat   fmri_data object with .dat (voxels x images) and .metadata_table
%         containing the subject_var column and all predictor/interaction
%         columns.
%
% Optional name-value pairs
% -------------------------
%   'predictors'         cellstr of metadata column names. Each becomes a
%                        binary 'Same<col>' column in the output table.
%   'interactions'       cell of cellstr (each inner cell length 2). Each
%                        pair-of-predictors becomes an element-wise AND
%                        interaction column.
%   'three_way'          cell of cellstr (each inner cell length 3). 3-way ANDs.
%   'subject_var'        Metadata column for subject IDs. Default 'subject_id'.
%   'pair_scope'         'within_subject' (default) | 'all'. The LME path
%                        wants within-subject. The fitlm/rsa_lm fixed-effects
%                        path may want all pairs.
%   'response_transform' 'fisherz' (default) | 'none' | 'rank'.
%   'metric'             RSM construction metric. Default 'correlation'.
%   'short_names'        struct mapping long metadata column -> short name
%                        for the output table (e.g. struct('subject_id','Subject')).
%                        Default: shortens commonly-suffixed names (drops _id,
%                        _number).
%   'verbose'            logical (default true).
%
% Outputs
% -------
%   tbl   table with columns:
%           Y                          double, similarity per pair
%           Subject                    categorical, grouping variable
%           Same<col>                  binary 0/1 per predictor
%           Same<col1>x<col2>          binary 0/1 per interaction
%
%   info  struct:
%           .predictor_names           cellstr of Same<col> column names
%           .interaction_names         cellstr
%           .three_way_names           cellstr
%           .subject_var               actual subject column name
%           .response_transform        applied
%           .n_pairs_per_subject       [n_subjects x 1]
%           .formula_skeleton          suggested Wilkinson formula stub
%
% References
% ----------
% Reproduces the per-subject table assembly in `08072024 Run-Level RDM
% Analysis with RSA Toolbox.mlx` lines 1804-1957.

% =========================================================================
% Parse inputs
% =========================================================================
p = inputParser;
p.addParameter('predictors',         {},             @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
p.addParameter('interactions',       {},             @iscell);
p.addParameter('three_way',          {},             @iscell);
p.addParameter('subject_var',        'subject_id',   @(x) ischar(x) || isstring(x));
p.addParameter('pair_scope',         'within_subject', @(x) (ischar(x) || isstring(x)) && ismember(lower(char(x)), {'within_subject','all'}));
p.addParameter('response_transform', 'fisherz',      @(x) (ischar(x) || isstring(x)) && ismember(lower(char(x)), {'fisherz','none','rank'}));
p.addParameter('metric',             'correlation',  @(x) (ischar(x) || isstring(x)) && ismember(lower(char(x)), {'correlation','spearman','cosine'}));
p.addParameter('short_names',        struct(),       @isstruct);
p.addParameter('verbose',            true,           @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

predictors   = cellstr(opt.predictors);
interactions = opt.interactions;
three_way    = opt.three_way;

% Validate
if isempty(predictors)
    error('assemble_lme_table:noPredictors', ...
        'Need at least one predictor metadata column.');
end

mt = dat.metadata_table;
if isempty(mt)
    error('assemble_lme_table:noMetadata', 'dat.metadata_table is empty.');
end

subject_var = char(opt.subject_var);
% Allow common subject column aliases
if ~ismember(subject_var, mt.Properties.VariableNames)
    aliases = {'subject_id','sub','subject','subj','sid'};
    hit = '';
    for i = 1:numel(aliases)
        if ismember(aliases{i}, mt.Properties.VariableNames), hit = aliases{i}; break; end
    end
    if isempty(hit)
        error('assemble_lme_table:noSubjectVar', ...
            ['No subject column found in metadata_table. Tried: %s. ', ...
             'Pass ''subject_var'' explicitly.'], strjoin([{subject_var}, aliases], ', '));
    end
    if opt.verbose
        fprintf('assemble_lme_table: subject_var=''%s'' not found; using ''%s''.\n', subject_var, hit);
    end
    subject_var = hit;
end

% Subject-in-predictors guard (design doc §6.2): only for the within-subject
% (LME) path, where Subject is the random-effects grouping. For pair_scope=
% 'all' (the rsa_lm fixed-effects path) SameSubject is a valid predictor.
if strcmpi(opt.pair_scope, 'within_subject') && any(strcmp(predictors, subject_var))
    error('assemble_lme_table:subjectInPredictors', ...
        ['"%s" is the random-effects grouping variable; remove it from ''predictors''. ', ...
         'Use rsa_lm() if you want Subject as a fixed effect instead.'], subject_var);
end

% For the all-pairs (fixed-effects) path, SameSubject is a valid predictor.
% Ensure subject_var is included exactly once in the predictor list so it is
% processed uniformly and reflected in info.predictor_names (no separate
% special-casing, no duplication).
if strcmpi(opt.pair_scope, 'all') && ~any(strcmp(predictors, subject_var))
    predictors = [{subject_var}, predictors];
end

% Check all predictor columns exist
for i = 1:numel(predictors)
    if ~ismember(predictors{i}, mt.Properties.VariableNames)
        error('assemble_lme_table:missingCol', ...
            'Predictor "%s" not found in metadata_table.', predictors{i});
    end
end

% Build short-name mapping (short version used as Y-column suffix)
short_names = opt.short_names;
for i = 1:numel(predictors)
    if ~isfield(short_names, predictors{i})
        short_names.(predictors{i}) = default_short_name(predictors{i});
    end
end
if ~isfield(short_names, subject_var)
    short_names.(subject_var) = default_short_name(subject_var);
end

% =========================================================================
% Per-subject table blocks
% =========================================================================
sub_ids_all = mt.(subject_var);
[~, subject_unique] = group_unique(sub_ids_all);
n_subj = numel(subject_unique);

block_tables = cell(n_subj, 1);
n_pairs_per_subject = zeros(n_subj, 1);

X = double(dat.dat);   % voxels x n_images

if strcmpi(opt.pair_scope, 'within_subject')
    % Within-subject pairs only (default LME path)
    for s = 1:n_subj
        block_tables{s} = build_subject_block( ...
            X, mt, sub_ids_all, subject_unique{s}, subject_var, ...
            predictors, interactions, three_way, short_names, opt);
        n_pairs_per_subject(s) = height(block_tables{s});
    end
else
    % All pairs (fitlm path) -- treat as one giant block
    block_tables{1} = build_omnibus_block( ...
        X, mt, sub_ids_all, subject_var, ...
        predictors, interactions, three_way, short_names, opt);
    n_pairs_per_subject = height(block_tables{1});
end

tbl = vertcat(block_tables{:});

% =========================================================================
% Build info struct + formula skeleton
% =========================================================================
info = struct();
info.predictor_names       = cellfun(@(c) ['Same' short_names.(c)], predictors, 'UniformOutput', false);
info.interaction_names     = cellfun(@(p) interaction_name(p, short_names), interactions, 'UniformOutput', false);
info.three_way_names       = cellfun(@(p) interaction_name(p, short_names), three_way,    'UniformOutput', false);
info.subject_var           = subject_var;
info.subject_var_short     = short_names.(subject_var);
info.response_transform    = opt.response_transform;
info.n_pairs_per_subject   = n_pairs_per_subject;
info.metric                = opt.metric;
info.pair_scope            = opt.pair_scope;

% Suggested formula skeleton
rhs_terms = [info.predictor_names, info.interaction_names, info.three_way_names];
if strcmpi(opt.pair_scope, 'within_subject')
    info.formula_skeleton = sprintf('Y ~ %s + (1 | %s)', ...
        strjoin(rhs_terms, ' + '), info.subject_var_short);
else
    info.formula_skeleton = sprintf('Y ~ %s', strjoin(rhs_terms, ' + '));
end

if opt.verbose
    fprintf('assemble_lme_table: %d rows, %d subjects, %d predictors, %d interactions\n', ...
        height(tbl), n_subj, numel(predictors), numel(interactions));
end

end


% =========================================================================
function block = build_subject_block(X, mt, sub_ids_all, sid, subject_var, ...
                                     predictors, interactions, three_way, short_names, opt)

is_s = match_subject(sub_ids_all, sid);
sub_X = X(:, is_s);
sub_meta = mt(is_s, :);

R_s = compute_rsm_single(sub_X, opt.metric);
% Upper-triangle vector
k = size(R_s, 1);
mask = triu(true(k), 1);
y = R_s(mask);

% Apply response transform
y = apply_response_transform(y, opt.response_transform);

n_pairs = numel(y);
block = table(y, 'VariableNames', {'Y'});

% Subject column (categorical)
sub_short = short_names.(subject_var);
block.(sub_short) = categorical(repmat({char(string(sid))}, n_pairs, 1));

% Predictor columns: SamePredictor = (v(i) == v(j)) on the upper-tri pairs
pred_vecs = struct();   % store for interactions
for p = 1:numel(predictors)
    col = predictors{p};
    v = sub_meta.(col);
    M = same_value_matrix(v);
    vec = M(mask);
    short_col = ['Same' short_names.(col)];
    block.(short_col) = double(vec);
    pred_vecs.(col) = double(vec);
end

% Interaction columns: element-wise AND of constituent predictors
for k_i = 1:numel(interactions)
    parts = interactions{k_i};
    name = interaction_name(parts, short_names);
    v = ones(n_pairs, 1);
    for p = 1:numel(parts)
        v = v .* pred_vecs.(parts{p});
    end
    block.(name) = double(v);
end

% Three-way interactions
for k_i = 1:numel(three_way)
    parts = three_way{k_i};
    name = interaction_name(parts, short_names);
    v = ones(n_pairs, 1);
    for p = 1:numel(parts)
        v = v .* pred_vecs.(parts{p});
    end
    block.(name) = double(v);
end

end


% =========================================================================
function block = build_omnibus_block(X, mt, sub_ids_all, subject_var, ...
                                     predictors, interactions, three_way, short_names, opt)
% Build a single all-pairs table (used for rsa_lm / pair_scope='all').

R_all = compute_rsm_single(X, opt.metric);
k = size(R_all, 1);
mask = triu(true(k), 1);
y = apply_response_transform(R_all(mask), opt.response_transform);
n_pairs = numel(y);

block = table(y, 'VariableNames', {'Y'});

% subject_var has already been folded into `predictors` by the caller, so it
% is built here as just another Same<X> predictor -- no special-casing.
pred_vecs = struct();
for p = 1:numel(predictors)
    col = predictors{p};
    v = mt.(col);
    M = same_value_matrix(v);
    vec = M(mask);
    short_col = ['Same' short_names.(col)];
    block.(short_col) = double(vec);
    pred_vecs.(col) = double(vec);
end

% Interactions
for k_i = 1:numel(interactions)
    parts = interactions{k_i};
    name = interaction_name(parts, short_names);
    v = ones(n_pairs, 1);
    for p = 1:numel(parts), v = v .* pred_vecs.(parts{p}); end
    block.(name) = double(v);
end
for k_i = 1:numel(three_way)
    parts = three_way{k_i};
    name = interaction_name(parts, short_names);
    v = ones(n_pairs, 1);
    for p = 1:numel(parts), v = v .* pred_vecs.(parts{p}); end
    block.(name) = double(v);
end

end


% =========================================================================
function R = compute_rsm_single(X, metric)
% Compute one similarity matrix from voxels x images data.
m = lower(char(metric));
switch m
    case 'correlation'
        R = corr(X, 'Type', 'Pearson', 'rows', 'pairwise');
    case 'spearman'
        R = corr(X, 'Type', 'Spearman', 'rows', 'pairwise');
    case 'cosine'
        norms = sqrt(sum(X.^2, 1));
        denom = norms' * norms;
        R = (X' * X) ./ denom;
    otherwise
        error('assemble_lme_table:badMetric', 'Unknown metric: %s', m);
end
R = (R + R') / 2;   % symmetrize
end


function y = apply_response_transform(y, mode)
mode = lower(char(mode));
switch mode
    case 'fisherz'
        y(y >  0.9999999) =  0.9999999;
        y(y < -0.9999999) = -0.9999999;
        y = atanh(y);
    case 'rank'
        y = tiedrank(y);
    case 'none'
        % no-op
end
end


function M = same_value_matrix(v)
% k x k binary matrix where (i, j) = 1 iff v(i) == v(j).
if iscell(v) || isstring(v) || iscategorical(v)
    [~, ~, codes] = unique(v, 'stable');
else
    codes = double(v);
end
codes = codes(:);
M = double(codes == codes');
end


function name = interaction_name(parts, short_names)
% Predictable interaction column name: join the individual Same<X> column
% names with 'x'. E.g. {condition,bodysite} -> 'SameConditionxSameBodysite'.
% This lets users construct the interaction column name directly from the
% main-effect column names they already know.
sames = cellfun(@(c) ['Same' short_names.(c)], parts, 'UniformOutput', false);
name = strjoin(sames, 'x');
end


function s = default_short_name(col)
% Drop common suffixes; capitalize first letter.
s = char(col);
s = regexprep(s, '_id$',     '');
s = regexprep(s, '_number$', '');
s = regexprep(s, '_no$',     '');
if isempty(s), s = char(col); end
if isstrprop(s(1), 'lower'), s(1) = upper(s(1)); end
% strip remaining underscores (won't be valid in Wilkinson formula without escapes)
s = regexprep(s, '[^A-Za-z0-9]', '');
end


function is_match = match_subject(v, sid)
sid = char(string(sid));
if iscell(v), is_match = strcmp(v, sid);
elseif iscategorical(v) || isstring(v), is_match = strcmp(string(v), sid);
else, is_match = (double(v) == double(string(sid)));
end
end


function [G, vals] = group_unique(v)
[G, vals] = findgroups(v);
if iscategorical(vals), vals = cellstr(vals);
elseif isstring(vals), vals = cellstr(vals);
elseif iscell(vals), % already cellstr
else, vals = cellfun(@(x) char(string(x)), num2cell(vals), 'UniformOutput', false);
end
end
