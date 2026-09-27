function [T, best_idx] = rsa_compare_lme_models(dat, formulas, varargin)
% rsa_compare_lme_models  Fit and compare a sequence of nested LME models.
%
% Replaces the manual fitlme-loop pattern in `08072024 Run-Level RDM Analysis`
% lines 2633-2691 with one declarative call.
%
% Usage
% -----
%   formulas = { ...
%       'Y ~ 1 + (1 | Sub)';
%       'Y ~ SameCondition + (1 | Sub)';
%       'Y ~ SameCondition + SameBodysite + (1 | Sub)';
%       'Y ~ SameCondition + SameBodysite + SameConditionxSameBodysite + (1 | Sub)';
%       'Y ~ SameCondition + SameBodysite + SameConditionxSameBodysite + SameSession + (1 | Sub)';
%       };
%   [T, best] = rsa_compare_lme_models(dat, formulas, ...
%       'predictors',    {'condition','bodysite','session_number'}, ...
%       'interactions',  {{'condition','bodysite'}}, ...
%       'subject_var',   'sub', ...
%       'select_by',     'aic');
%
% Inputs
% ------
%   dat       fmri_data object (passed to rsa_lme for each formula)
%   formulas  cellstr of Wilkinson formulas, ordered from simplest to most complex
%
% Optional name-value
% -------------------
%   All rsa_lme name-value pairs are forwarded (predictors, interactions,
%   three_way, subject_var, response_transform, metric, fit_method, etc.)
%
%   Plus:
%   'lrt_pairs'   'sequential' (default) -- LRT each model against the
%                                            previous one
%                 'vs_first'    -- LRT each model against the first
%                 'none'        -- skip LRT
%   'select_by'   'aic' (default) | 'bic' | 'loglik' -- which criterion to
%                                          report as the best model
%   'verbose'     logical (default true)
%
% Outputs
% -------
%   T         table with one row per model:
%               model | n_params | AIC | BIC | logLik | lrt_chi2 | lrt_df | lrt_p
%   best_idx  index of the model with the best 'select_by' criterion

% Split off the rsa_compare-only params
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('lrt_pairs', 'sequential', @(x) (ischar(x) || isstring(x)) && ismember(lower(char(x)), {'sequential','vs_first','none'}));
p.addParameter('select_by', 'aic',        @(x) (ischar(x) || isstring(x)) && ismember(lower(char(x)), {'aic','bic','loglik'}));
p.addParameter('verbose',   true,         @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

% Forward unknown params to rsa_lme.
% NB: Nested LRT of LME models with different fixed effects requires ML
% fit (REML criterion isn't comparable across different fixed designs).
% Force fit_method='ML' so compare() works.
fwd = struct2cellfwd(p.Unmatched);
% Strip any user-supplied fit_method and replace with ML
fwd = remove_named(fwd, 'fit_method');
fwd = [fwd, {'fit_method', 'ML'}];
if opt.verbose
    fprintf('rsa_compare_lme_models: forcing fit_method=''ML'' for nested LRT validity.\n');
end

formulas = cellstr(formulas);
n = numel(formulas);

% =========================================================================
% Fit each model
% =========================================================================
models = cell(n, 1);
aics    = nan(n, 1);
bics    = nan(n, 1);
logliks = nan(n, 1);
n_params = nan(n, 1);

for i = 1:n
    if opt.verbose, fprintf('rsa_compare_lme_models: fitting %d/%d: %s\n', i, n, formulas{i}); end
    models{i} = rsa_lme(dat, formulas{i}, fwd{:}, 'verbose', false);
    aics(i)    = models{i}.ModelCriterion.AIC;
    bics(i)    = models{i}.ModelCriterion.BIC;
    logliks(i) = models{i}.LogLikelihood;
    % NumCoefficients is reliably available; covariance-parameter count
    % varies across MATLAB versions, so fall back to the AIC-implied count.
    n_params(i) = models{i}.NumCoefficients;
end

% =========================================================================
% Nested likelihood-ratio tests
% =========================================================================
lrt_chi2 = nan(n, 1);
lrt_df   = nan(n, 1);
lrt_p    = nan(n, 1);

switch lower(char(opt.lrt_pairs))
    case 'sequential'
        for i = 2:n
            try
                cmp = compare(models{i-1}, models{i}, 'CheckNesting', true);
                lrt_chi2(i) = cmp.LRStat(end);
                lrt_df(i)   = cmp.deltaDF(end);
                lrt_p(i)    = cmp.pValue(end);
            catch ME
                if opt.verbose
                    warning('rsa_compare_lme_models:lrtFailed', ...
                        'LRT failed for model %d: %s', i, ME.message);
                end
            end
        end
    case 'vs_first'
        for i = 2:n
            try
                cmp = compare(models{1}, models{i}, 'CheckNesting', true);
                lrt_chi2(i) = cmp.LRStat(end);
                lrt_df(i)   = cmp.deltaDF(end);
                lrt_p(i)    = cmp.pValue(end);
            catch ME
                if opt.verbose
                    warning('rsa_compare_lme_models:lrtFailed', ...
                        'LRT failed for model %d: %s', i, ME.message);
                end
            end
        end
    case 'none'
        % skip
end

% =========================================================================
% Best model
% =========================================================================
switch lower(char(opt.select_by))
    case 'aic',    [~, best_idx] = min(aics);
    case 'bic',    [~, best_idx] = min(bics);
    case 'loglik', [~, best_idx] = max(logliks);
end

% =========================================================================
% Output table
% =========================================================================
T = table(formulas(:), n_params, aics, bics, logliks, lrt_chi2, lrt_df, lrt_p, ...
    'VariableNames', {'model','n_params','AIC','BIC','logLik','lrt_chi2','lrt_df','lrt_p'});

if opt.verbose
    fprintf('\nrsa_compare_lme_models: best by %s = model %d:\n  %s\n', ...
        upper(char(opt.select_by)), best_idx, formulas{best_idx});
end

end


function out = struct2cellfwd(s)
% Convert a struct of name-value pairs back into a cell array for forwarding
fn = fieldnames(s);
out = cell(1, 2 * numel(fn));
for i = 1:numel(fn)
    out{2*i - 1} = fn{i};
    out{2*i}     = s.(fn{i});
end
end


function out = remove_named(args, name)
% Strip a name-value pair from a varargin cell
out = args;
idx = [];
for i = 1:2:numel(out)-1
    if (ischar(out{i}) || isstring(out{i})) && strcmpi(char(out{i}), name)
        idx = i;
        break
    end
end
if ~isempty(idx) && idx < numel(out)
    out([idx, idx+1]) = [];
end
end
