function [T, mdls, info] = rsa_lm_by_subject(dat, varargin)
% rsa_lm_by_subject  Per-subject fitlm fits with aggregated coefficient table.
%
% Replicates the per-subject loop from `08072024 Run-Level RDM Analysis with
% RSA Toolbox.mlx` lines 1879-1900 + lines 2050-2079 (partial R²): for each
% subject, fits a separate `fitlm` of similarity ~ same-vs-different predictors
% on that subject's within-subject pairs only. Returns a long-format table of
% coefficients with subject IDs so you can run group-level inference on the
% per-subject betas (paired ttest across subjects, etc.).
%
% Differs from rsa_lme: per-subject FIXED-effects fits aggregated post hoc,
% rather than one mixed-effects model treating Subject as a random effect.
% Useful for assessing between-subject variability in coefficients and for
% sanity-checking the LME estimates.
%
% Usage
% -----
%   T = rsa_lm_by_subject(dat, ...
%       'predictors',   {'condition','bodysite','sesno'}, ...
%       'interactions', {{'condition','bodysite'}}, ...
%       'subject_var',  'sub');
%
% Optional name-value
% -------------------
%   Same as rsa_lme except no random-effects spec. Plus:
%   'partial_r2'    logical (default true). Compute per-term partial R² by
%                   refitting the reduced model dropping each term. Adds
%                   columns to T.
%   'verbose'       logical (default true).
%
% Outputs
% -------
%   T     long-format table:
%           sub | term | beta | se | t | p | partial_R2 | full_R2 | full_adj_R2
%   mdls  cell array of fitted LinearModel objects, one per subject
%   info  struct from the underlying assemble_lme_table call (uses the
%         first subject's metadata as representative)

p = inputParser;
p.addParameter('predictors',         {},             @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
p.addParameter('interactions',       {},             @iscell);
p.addParameter('three_way',          {},             @iscell);
p.addParameter('subject_var',        'subject_id',   @(x) ischar(x) || isstring(x));
p.addParameter('response_transform', 'fisherz',      @(x) ischar(x) || isstring(x));
p.addParameter('metric',             'correlation',  @(x) ischar(x) || isstring(x));
p.addParameter('partial_r2',         true,           @(x) islogical(x) || isnumeric(x));
p.addParameter('verbose',            true,           @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

predictors = cellstr(opt.predictors);
if isempty(predictors)
    error('rsa_lm_by_subject:noPredictors', 'Pass at least one predictor.');
end

% Build the full within-subject pairs table; we'll slice by Subject below
[tbl_all, info] = assemble_lme_table(dat, ...
    'predictors',         predictors, ...
    'interactions',       opt.interactions, ...
    'three_way',          opt.three_way, ...
    'subject_var',        opt.subject_var, ...
    'pair_scope',         'within_subject', ...
    'response_transform', opt.response_transform, ...
    'metric',             opt.metric, ...
    'verbose',            opt.verbose);

sub_short = info.subject_var_short;
sub_levels = categories(tbl_all.(sub_short));
n_subj = numel(sub_levels);

% Predictor + interaction column names in the table
term_cols = [info.predictor_names, info.interaction_names, info.three_way_names];
n_terms = numel(term_cols);

% Build base formula (no random effects)
rhs = strjoin(term_cols, ' + ');
formula = sprintf('Y ~ %s', rhs);

% =========================================================================
% Per-subject fits
% =========================================================================
mdls = cell(n_subj, 1);
rows = cell(n_subj, 1);

for s = 1:n_subj
    is_s = tbl_all.(sub_short) == sub_levels{s};
    tbl_s = tbl_all(is_s, :);
    mdls{s} = fitlm(tbl_s, formula);

    n_coefs = mdls{s}.NumCoefficients;
    sub_col   = repmat({char(sub_levels{s})}, n_coefs, 1);
    coef_name = mdls{s}.Coefficients.Properties.RowNames;
    betas     = mdls{s}.Coefficients.Estimate;
    ses       = mdls{s}.Coefficients.SE;
    ts        = mdls{s}.Coefficients.tStat;
    ps        = mdls{s}.Coefficients.pValue;
    full_r2     = repmat(mdls{s}.Rsquared.Ordinary, n_coefs, 1);
    full_adj_r2 = repmat(mdls{s}.Rsquared.Adjusted, n_coefs, 1);

    % Partial R² per term: refit dropping that term, R²_full - R²_reduced
    partial_r2 = nan(n_coefs, 1);
    if logical(opt.partial_r2)
        for c = 1:n_coefs
            name = coef_name{c};
            if strcmp(name, '(Intercept)'), continue; end
            others = setdiff(term_cols, {name});
            if isempty(others)
                f_red = 'Y ~ 1';
            else
                f_red = sprintf('Y ~ %s', strjoin(others, ' + '));
            end
            try
                mdl_red = fitlm(tbl_s, f_red);
                partial_r2(c) = mdls{s}.Rsquared.Ordinary - mdl_red.Rsquared.Ordinary;
            catch
                partial_r2(c) = NaN;
            end
        end
    end

    rows{s} = table(sub_col, coef_name, betas, ses, ts, ps, partial_r2, ...
        full_r2, full_adj_r2, ...
        'VariableNames', {sub_short, 'term', 'beta', 'se', 't', 'p', ...
                          'partial_R2', 'full_R2', 'full_adj_R2'});
end

T = vertcat(rows{:});

if opt.verbose
    fprintf('rsa_lm_by_subject: fit %d per-subject models with %d terms each.\n', ...
        n_subj, n_terms);
end

end
