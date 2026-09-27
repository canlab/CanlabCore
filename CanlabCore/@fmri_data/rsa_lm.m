function [mdl, tbl, info] = rsa_lm(dat, varargin)
% rsa_lm  Fixed-effects multi-level RSA via fitlm.
%
% Pools all (i, j) upper-triangle pairs from the omnibus image-level RSM
% (including cross-subject pairs by default per design doc §6.6) and fits
% a linear model with same-vs-different predictors. Subject is treated as
% a same-vs-different predictor itself (SameSubject), NOT as a random effect.
%
% Use rsa_lme() instead if you want random effects.
%
% Usage examples
% --------------
%   % Default: all pairs, all predictors
%   mdl = rsa_lm(dat, ...
%       'predictors',  {'subject_id','session_number','condition','bodysite'});
%
%   % Within-subject only (subset of the data; conceptually equivalent to
%   % the LME fixed-effect part)
%   mdl = rsa_lm(dat, ...
%       'predictors',  {'condition','bodysite','session_number'}, ...
%       'pair_scope',  'within_subject');
%
%   % With interactions
%   mdl = rsa_lm(dat, ...
%       'predictors',   {'condition','bodysite','session_number'}, ...
%       'interactions', {{'condition','bodysite'},{'session_number','condition'}});
%
% Optional name-value
% -------------------
%   'predictors'         cellstr of metadata columns.
%   'interactions'       cell of pair-cellstr.
%   'three_way'          cell of triple-cellstr.
%   'subject_var'        Default 'subject_id'.
%   'pair_scope'         'all' (default) | 'within_subject'.
%   'response_transform' 'fisherz' (default) | 'none' | 'rank'.
%   'metric'             'correlation' (default) | 'spearman' | 'cosine'.
%   'standardize'        logical (default false). Standardize predictors
%                        before fit.
%   'verbose'            logical (default true).
%
% Outputs
% -------
%   mdl   LinearModel object
%   tbl   long-format table used for fitting
%   info  struct from assemble_lme_table

p = inputParser;
p.addParameter('predictors',         {},             @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
p.addParameter('interactions',       {},             @iscell);
p.addParameter('three_way',          {},             @iscell);
p.addParameter('subject_var',        'subject_id',   @(x) ischar(x) || isstring(x));
p.addParameter('pair_scope',         'all',          @(x) (ischar(x) || isstring(x)) && ismember(lower(char(x)), {'all','within_subject'}));
p.addParameter('response_transform', 'fisherz',      @(x) ischar(x) || isstring(x));
p.addParameter('metric',             'correlation',  @(x) ischar(x) || isstring(x));
p.addParameter('standardize',        false,          @(x) islogical(x) || isnumeric(x));
p.addParameter('verbose',            true,           @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

predictors = cellstr(opt.predictors);
if isempty(predictors)
    error('rsa_lm:noPredictors', 'Pass at least one predictor.');
end

% =========================================================================
% Assemble table
% =========================================================================
[tbl, info] = assemble_lme_table(dat, ...
    'predictors',         predictors, ...
    'interactions',       opt.interactions, ...
    'three_way',          opt.three_way, ...
    'subject_var',        opt.subject_var, ...
    'pair_scope',         opt.pair_scope, ...
    'response_transform', opt.response_transform, ...
    'metric',             opt.metric, ...
    'verbose',            opt.verbose);

% =========================================================================
% Build formula
% =========================================================================
% All predictor + interaction Same<X> columns are the fixed effects.
% For pair_scope='all', subject_var was folded into predictors by
% assemble_lme_table, so SameSubject is already in info.predictor_names.
rhs_terms = [info.predictor_names, info.interaction_names, info.three_way_names];
rhs_terms = unique(rhs_terms, 'stable');   % guard against any duplication

% Optional standardize: z-score predictor columns in-place
if logical(opt.standardize)
    for i = 1:numel(rhs_terms)
        col = rhs_terms{i};
        v = tbl.(col);
        if std(v) > 0, tbl.(col) = (v - mean(v)) ./ std(v); end
    end
end

formula = sprintf('Y ~ %s', strjoin(rhs_terms, ' + '));

if opt.verbose
    fprintf('rsa_lm: fitting formula:\n  %s\n  (n=%d rows, scope=%s)\n', ...
        formula, height(tbl), opt.pair_scope);
end

% =========================================================================
% Fit
% =========================================================================
mdl = fitlm(tbl, formula);

if opt.verbose
    fprintf('rsa_lm: R^2 = %.4f, adjusted R^2 = %.4f\n', ...
        mdl.Rsquared.Ordinary, mdl.Rsquared.Adjusted);
end

end
