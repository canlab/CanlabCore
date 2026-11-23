function T = rsa_partial_r2(mdl, tbl)
% rsa_partial_r2  Partial R² for each predictor in a fitted LinearModel.
%
% For each non-intercept predictor, refits the model dropping that predictor
% and computes partial R² = full_R² - reduced_R². This is the standard
% "incremental variance" metric. Used in `08072024 Run-Level RDM Analysis`
% lines 2050-2079.
%
% Usage
% -----
%   [mdl_lm, tbl, ~] = dat.rsa_lm('predictors', {...});
%   T = rsa_partial_r2(mdl_lm, tbl);
%   % T columns: term | partial_R2 | Cohens_d
%
% Inputs
% ------
%   mdl  LinearModel object (from fitlm)
%   tbl  table used to fit mdl
%
% Output
% ------
%   T  table:
%       term        char
%       partial_R2  double
%       Cohens_d    double (per-predictor effect size on Y, computed from std)

assert(isa(mdl, 'LinearModel'), 'rsa_partial_r2: requires a LinearModel input.');
assert(istable(tbl), 'rsa_partial_r2: requires a table input.');

terms = mdl.Coefficients.Properties.RowNames;
n_terms = numel(terms);

partial_r2 = nan(n_terms, 1);
cohens_d   = nan(n_terms, 1);

% Identify the response column from the model formula
y_name = mdl.ResponseName;
std_y  = std(tbl.(y_name), 'omitnan');

% All predictor names from the model (excluding intercept)
all_predictors = setdiff(terms, {'(Intercept)'}, 'stable');

for i = 1:n_terms
    name = terms{i};
    if strcmp(name, '(Intercept)'), continue; end

    others = setdiff(all_predictors, {name}, 'stable');
    if isempty(others)
        f_red = sprintf('%s ~ 1', y_name);
    else
        f_red = sprintf('%s ~ %s', y_name, strjoin(others, ' + '));
    end
    try
        mdl_red = fitlm(tbl, f_red);
        partial_r2(i) = mdl.Rsquared.Ordinary - mdl_red.Rsquared.Ordinary;
    catch
        partial_r2(i) = NaN;
    end

    % Cohen's d on this predictor (beta normalized by Y std)
    if std_y > 0
        cohens_d(i) = mdl.Coefficients.Estimate(i) / std_y;
    end
end

T = table(terms, partial_r2, cohens_d, ...
    'VariableNames', {'term', 'partial_R2', 'Cohens_d'});

end
