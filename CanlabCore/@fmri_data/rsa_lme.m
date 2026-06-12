function [mdl, tbl, info] = rsa_lme(dat, varargin)
% rsa_lme  Random-effects multi-level RSA via fitlme.
%
% Implements the Phase 3 LME pipeline: assemble within-subject (i, j)
% upper-triangle pairs of the per-subject image-level RSM into a long-format
% table, then fit fitlme with the requested fixed and random effects.
%
% Two ways to specify the model
% -----------------------------
%   Form 1 -- Wilkinson formula (explicit):
%     mdl = rsa_lme(dat, 'Y ~ SameCondition + SameBodysite + (1 | Subject)', ...
%         'predictors', {'condition','bodysite'}, 'subject_var','sub');
%
%   Form 2 -- Name-value assembly (auto formula):
%     mdl = rsa_lme(dat, ...
%         'predictors',    {'condition','bodysite','session_number'}, ...
%         'fixed_effects', {'condition','bodysite','condition:bodysite'}, ...
%         'random_effects',{'(1 | sub)'}, ...
%         'subject_var',   'sub');
%
% Usage examples
% --------------
%   % Simple random-intercept model
%   mdl = rsa_lme(dat, ...
%       'predictors',    {'condition','bodysite','session_number'}, ...
%       'subject_var',   'sub');
%   % Auto formula: 'Y ~ SameCondition + SameBodysite + SameSession + (1 | Sub)'
%
%   % Add an interaction
%   mdl = rsa_lme(dat, ...
%       'predictors',    {'condition','bodysite'}, ...
%       'interactions',  {{'condition','bodysite'}}, ...
%       'subject_var',   'sub');
%
%   % Random slopes for SameCondition by subject
%   mdl = rsa_lme(dat, ...
%       'Y ~ SameCondition + SameBodysite + (SameCondition | Sub)', ...
%       'predictors', {'condition','bodysite'}, ...
%       'subject_var','sub');
%
% Optional name-value
% -------------------
%   'predictors'         cellstr of metadata columns to include as same-vs-
%                        different fixed-effect RDMs.
%   'fixed_effects'      cellstr of Wilkinson terms to include in the
%                        fixed-effect formula (used only with name-value
%                        form). Use ':' for interactions, e.g. 'condition:bodysite'.
%   'random_effects'     cellstr of random-effect Wilkinson terms, e.g.
%                        {'(1 | sub)'}, {'(condition | sub)'}.
%   'interactions'       cell of cellstr pairs to ALSO add as fixed-effect
%                        columns (built as element-wise AND).
%   'subject_var'        Default 'subject_id' (also matches 'sub', etc.).
%   'response_transform' 'fisherz' (default) | 'none' | 'rank'.
%   'fit_method'         'REML' (default) | 'ML'.
%   'verbose'            logical (default true).
%
% Outputs
% -------
%   mdl   LinearMixedModel object
%   tbl   the long-format table used for fitting (height = sum of
%         within-subject pair counts)
%   info  struct from assemble_lme_table

% =========================================================================
% Detect form 1 (formula) vs form 2 (name-value)
% =========================================================================
formula = '';
if ~isempty(varargin) && (ischar(varargin{1}) || isstring(varargin{1}))
    first = char(varargin{1});
    % Treat as formula if it contains '~' AND isn't a known parameter name
    if contains(first, '~')
        formula = first;
        varargin = varargin(2:end);
    end
end

% =========================================================================
% Parse name-value args
% =========================================================================
p = inputParser;
p.addParameter('predictors',         {},             @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
p.addParameter('fixed_effects',      {},             @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
p.addParameter('random_effects',     {},             @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
p.addParameter('interactions',       {},             @iscell);
p.addParameter('three_way',          {},             @iscell);
p.addParameter('subject_var',        'subject_id',   @(x) ischar(x) || isstring(x));
p.addParameter('response_transform', 'fisherz',      @(x) ischar(x) || isstring(x));
p.addParameter('metric',             'correlation',  @(x) ischar(x) || isstring(x));
p.addParameter('fit_method',         'REML',         @(x) (ischar(x) || isstring(x)) && ismember(upper(char(x)), {'REML','ML'}));
p.addParameter('verbose',            true,           @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

predictors   = cellstr(opt.predictors);
fixed_eff    = cellstr(opt.fixed_effects);
random_eff   = cellstr(opt.random_effects);
interactions = opt.interactions;
three_way    = opt.three_way;

% =========================================================================
% Assemble table
% =========================================================================
[tbl, info] = assemble_lme_table(dat, ...
    'predictors',         predictors, ...
    'interactions',       interactions, ...
    'three_way',          three_way, ...
    'subject_var',        opt.subject_var, ...
    'pair_scope',         'within_subject', ...
    'response_transform', opt.response_transform, ...
    'metric',             opt.metric, ...
    'verbose',            opt.verbose);

% =========================================================================
% Build / validate formula
% =========================================================================
if isempty(formula)
    formula = build_formula_from_namevalue(info, predictors, fixed_eff, ...
        random_eff, interactions, three_way);
end

% Normalize + validate formula column names against the assembled table.
% Resolves both interaction-naming forms (SameAxSameB <-> SameAxB) and gives
% a clear error listing available columns if a term can't be matched.
formula = normalize_and_validate_formula(formula, tbl, info);

if opt.verbose
    fprintf('rsa_lme: fitting formula:\n  %s\n  (n=%d rows, %d subjects)\n', ...
        formula, height(tbl), numel(unique(tbl.(info.subject_var_short))));
end

% =========================================================================
% Fit
% =========================================================================
mdl = fitlme(tbl, formula, 'FitMethod', upper(char(opt.fit_method)));

if opt.verbose
    fprintf('rsa_lme: LME fit complete. AIC=%.1f, BIC=%.1f, logLik=%.1f\n', ...
        mdl.ModelCriterion.AIC, mdl.ModelCriterion.BIC, mdl.LogLikelihood);
end

end


function formula = build_formula_from_namevalue(info, predictors, fixed_eff, ...
                                                 random_eff, interactions, three_way) %#ok<INUSL>
% Convert name-value spec to a Wilkinson formula. Each predictor name is
% translated to its 'Same<short>' column name automatically.

% Map original metadata column name -> Same<Short> table column name
name_map = containers.Map(predictors, info.predictor_names);

% Determine fixed-effect terms
if isempty(fixed_eff)
    fixed_rhs = [info.predictor_names, info.interaction_names, info.three_way_names];
else
    fixed_rhs = cellfun(@(t) translate_term(t, name_map), fixed_eff, 'UniformOutput', false);
end
rhs_fixed = strjoin(fixed_rhs, ' + ');

% Random effects: default to (1 | Subject) if none specified
if isempty(random_eff)
    rhs_rand = sprintf('(1 | %s)', info.subject_var_short);
else
    rhs_rand = strjoin(cellfun(@(r) translate_random_effect(r, name_map, info.subject_var_short), ...
        random_eff, 'UniformOutput', false), ' + ');
end

formula = sprintf('Y ~ %s + %s', rhs_fixed, rhs_rand);
end


function out = translate_term(term, name_map)
% Translate a user-friendly term like 'condition' or 'condition:bodysite'
% to its same-vs-different column name 'SameCondition' / 'SameConditionxSameBodysite'.
term = strtrim(char(term));
if contains(term, ':')
    parts = strsplit(term, ':');
    parts = strtrim(parts);
    sames = cellfun(@(p) get_same_name(p, name_map), parts, 'UniformOutput', false);
    % Interaction column name = join the Same<X> names with 'x'
    out = strjoin(sames, 'x');
elseif contains(term, '*')
    parts = strsplit(term, '*');
    parts = strtrim(parts);
    % Main effects + 2-way interaction (Wilkinson 'a*b' = 'a + b + a:b')
    sames = cellfun(@(p) get_same_name(p, name_map), parts, 'UniformOutput', false);
    inter = strjoin(sames, 'x');
    out = strjoin([sames, {inter}], ' + ');
else
    out = get_same_name(term, name_map);
end
end


function out = get_same_name(name, name_map)
if isKey(name_map, name), out = name_map(name);
elseif startsWith(name, 'Same'), out = name;
else, out = ['Same' name];   % already a short form
end
end


function formula = normalize_and_validate_formula(formula, tbl, info)
% Check each fixed-effect term in the formula against the actual table
% columns. Resolve interaction-naming variants (SameAxSameB <-> SameAxB).
% If a term still can't be matched, error with the list of available columns.

actual_cols = tbl.Properties.VariableNames;

% Split formula into LHS ~ RHS
tilde = strfind(formula, '~');
if isempty(tilde)
    return   % not a formula we can parse; let fitlme handle it
end
lhs = strtrim(formula(1:tilde(1)-1));
rhs = strtrim(formula(tilde(1)+1:end));

% Pull out random-effect parenthetical clauses; leave them untouched
rand_clauses = regexp(rhs, '\([^)]*\)', 'match');
fixed_part   = regexprep(rhs, '\([^)]*\)', '');

% Split fixed terms on '+'
raw_terms = strtrim(strsplit(fixed_part, '+'));
raw_terms = raw_terms(~cellfun('isempty', raw_terms));

resolved   = {};
unresolved = {};
for i = 1:numel(raw_terms)
    t = raw_terms{i};
    if strcmp(t, '1') || ~isempty(regexp(t, '^[0-9.]+$', 'once'))
        resolved{end+1} = t; %#ok<AGROW>
        continue
    end
    col = resolve_column(t, actual_cols);
    if isempty(col)
        unresolved{end+1} = t; %#ok<AGROW>
    else
        resolved{end+1} = col; %#ok<AGROW>
    end
end

if ~isempty(unresolved)
    % Build a helpful error listing the predictor columns that DO exist
    pred_cols = [info.predictor_names, info.interaction_names, info.three_way_names];
    error('rsa_lme:unknownFormulaTerm', ...
        ['Formula term(s) not found in the assembled table: %s\n', ...
         'Available predictor columns are:\n  %s\n', ...
         'Interaction columns join the main-effect names with ''x'', ', ...
         'e.g. SameConditionxSameBodysite.\n', ...
         'Did you list all needed columns in ''predictors'' and ''interactions''?'], ...
        strjoin(unresolved, ', '), strjoin(pred_cols, ', '));
end

% Reassemble
rhs_new = strjoin(resolved, ' + ');
if ~isempty(rand_clauses)
    rhs_new = [rhs_new ' + ' strjoin(rand_clauses, ' + ')];
end
formula = sprintf('%s ~ %s', lhs, rhs_new);
end


function col = resolve_column(token, actual_cols)
% Return the actual column name matching `token`, trying interaction-naming
% variants. Returns '' if no match.
token = strtrim(token);
if ismember(token, actual_cols), col = token; return; end
if contains(token, 'x')
    % Variant A: collapse redundant 'Same' after each 'x' (SameAxSameB -> SameAxB)
    collapsed = regexprep(token, 'xSame', 'x');
    if ismember(collapsed, actual_cols), col = collapsed; return; end
    % Variant B: add 'Same' after each 'x' that lacks it (SameAxB -> SameAxSameB)
    expanded = regexprep(token, 'x(?!Same)', 'xSame');
    if ismember(expanded, actual_cols), col = expanded; return; end
end
col = '';
end


function out = translate_random_effect(re, name_map, subject_short)
% Random-effect term: '(1 | sub)' or '(condition | sub)' -- translate the
% grouping variable to its short name and any inner predictor to Same<X>.
re = char(re);
% Replace 'condition' with 'SameCondition' etc. inside the parens
keys = name_map.keys;
for i = 1:numel(keys)
    re = regexprep(re, ['(?<![a-zA-Z_])', keys{i}, '(?![a-zA-Z_])'], name_map(keys{i}));
end
% Translate subject grouping variable (right of |) -- common aliases
re = regexprep(re, '\|\s*subject_id\s*\)',  ['| ' subject_short ')']);
re = regexprep(re, '\|\s*sub\s*\)',         ['| ' subject_short ')']);
re = regexprep(re, '\|\s*subject\s*\)',     ['| ' subject_short ')']);
out = re;
end
