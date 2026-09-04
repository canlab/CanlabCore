function icc_struct = rsa_lme_icc(mdl)
% rsa_lme_icc  Extract variance-components ICC(s) from a fitted LinearMixedModel.
%
% For an LME of the form Y ~ ... + (1 | Group) + (X | Group) + ..., this
% returns the proportion of total variance explained by each random-effect
% variance component:
%
%   ICC_intercept = sigma_intercept^2 / total_variance
%   ICC_slope_X   = sigma_slope_X^2  / total_variance
%
% where total_variance = sum(random-effect variances) + residual variance.
%
% Reproduces the per-source variance calculation in `08072024 Run-Level RDM
% Analysis with RSA Toolbox.mlx` lines 2186-2225.
%
% Usage
% -----
%   mdl = dat.rsa_lme('Y ~ SameCondition + (SameCondition | Sub)', ...);
%   icc = rsa_lme_icc(mdl);
%   icc.summary             % table of source | variance | proportion (= ICC)
%   icc.total_variance      % scalar
%   icc.residual_variance   % scalar
%
% Output
% ------
%   icc_struct  struct with fields:
%     .summary           table: Source | Variance | ICC
%     .total_variance    sum of all variance components
%     .residual_variance scalar
%     .grouping_var      name of the random-effects grouping variable

assert(isa(mdl, 'LinearMixedModel'), 'rsa_lme_icc: requires a LinearMixedModel input.');

% Extract variance components
[psi, mse, stats] = covarianceParameters(mdl);
% psi is cell array of covariance matrices, one per random-effect group
% mse is the residual variance (scalar)
% stats has detail per group

% Build per-source list. psi{g}(i, i) is the variance of the i-th random
% effect in group g. We pull the term names from stats{g} where Type='std'
% and Name1==Name2 (diagonal entries).
sources    = {};
variances  = [];

for g = 1:numel(psi)
    grp_stats = stats{g};
    % Get group name (e.g. 'Sub'). Different MATLAB versions expose Group
    % as char matrix (one row per RE term), cellstr, or categorical.
    grp_name = sprintf('Group%d', g);
    try
        gv = grp_stats.Group;
        if ischar(gv)
            if size(gv, 1) >= 1, grp_name = strtrim(gv(1, :));   % char matrix
            else, grp_name = gv;
            end
        elseif iscell(gv) && ~isempty(gv), grp_name = char(gv{1});
        elseif iscategorical(gv) && ~isempty(gv), grp_name = char(gv(1));
        elseif isstring(gv) && ~isempty(gv), grp_name = char(gv(1));
        end
    catch
    end

    % Pull RE term names from Type='std' rows where Name1==Name2. Handle
    % both char-matrix and cellstr forms.
    re_names = arrayfun(@(i) sprintf('RE%d', i), 1:size(psi{g}, 1), 'UniformOutput', false);
    try
        types  = coerce_to_cellstr(grp_stats.Type);
        names1 = coerce_to_cellstr(grp_stats.Name1);
        names2 = coerce_to_cellstr(grp_stats.Name2);
        is_std  = strcmp(types, 'std');
        is_diag = strcmp(names1, names2);
        diag_rows = find(is_std & is_diag);
        if ~isempty(diag_rows)
            re_names = names1(diag_rows);
        end
    catch
        % keep generic fallback
    end

    n_re_terms = size(psi{g}, 1);
    for r = 1:n_re_terms
        if r <= numel(re_names), term_name = re_names{r};
        else, term_name = sprintf('RE%d', r); end
        sources{end+1, 1} = sprintf('%s | %s', term_name, grp_name); %#ok<AGROW>
        variances(end+1, 1) = psi{g}(r, r); %#ok<AGROW>
    end
end

% Add residual
sources{end+1, 1}   = 'Residual';
variances(end+1, 1) = mse;

total = sum(variances);
iccs  = variances ./ total;

icc_struct = struct();
icc_struct.summary           = table(sources, variances, iccs, ...
    'VariableNames', {'Source', 'Variance', 'ICC'});
icc_struct.total_variance    = total;
icc_struct.residual_variance = mse;

% Grouping variable(s) for convenience
grp_set = {};
for i = 1:numel(sources)
    tok = regexp(sources{i}, '\|\s*(\w+)', 'tokens', 'once');
    if ~isempty(tok), grp_set{end+1} = tok{1}; end %#ok<AGROW>
end
icc_struct.grouping_var = unique(grp_set);

end


function c = coerce_to_cellstr(v)
% Coerce a char matrix / cellstr / categorical / string to a cellstr column
if iscell(v), c = v(:);
elseif ischar(v)
    % rows of a char matrix become cells
    c = cell(size(v, 1), 1);
    for r = 1:size(v, 1), c{r} = strtrim(v(r, :)); end
elseif iscategorical(v), c = cellstr(v);
elseif isstring(v), c = cellstr(v);
else, c = cellstr(string(v));
end

end
