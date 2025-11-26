function [h, pairwise_tbl] = pairwise_contrasts(lme, varname, sig)
%PAIRWISE_CONTRASTS Pairwise tests for a categorical fixed effect in a LinearMixedModel.
%
% [h, pairwise_tbl] = pairwise_contrasts(lme, varname, sig)
%
% lme      : fitted LinearMixedModel (from fitlme)
% varname  : name of the categorical variable to compare (e.g., 'body_site')
% sig      : significance / correction option:
%              'none'        -> uncorrected p-values, alpha = 0.05
%              'fdr'         -> FDR-corrected p-values (Bh-FDR), alpha = 0.05 (default)
%              'bon'         -> Bonferroni-corrected p-values, alpha = 0.05
%              numeric alpha -> custom alpha (FDR by default)
%
% Returns:
%   h            : figure handle for the heatmap
%   pairwise_tbl : table of pairwise comparisons with F, df, p, p_fdr, p_bon
%
% Assumes default reference (dummy) coding for categorical variables in fitlme.

% Michael Sun, Ph.D. 11/26/2025 (Happy Thanksgiving-Eve)

    % ---------- Input handling ----------
    if nargin < 2 || isempty(varname)
        error('You must provide the name of a categorical variable (e.g., ''body_site'').');
    end

    if nargin < 3 || isempty(sig)
        sig_method = 'fdr';
        alpha = 0.05;
    elseif ischar(sig) || isstring(sig)
        sig = lower(string(sig));
        switch sig
            case "none"
                sig_method = 'none';
                alpha = 0.05;
            case "fdr"
                sig_method = 'fdr';
                alpha = 0.05;
            case {"bon","bonf","bonferroni"}
                sig_method = 'bon';
                alpha = 0.05;
            otherwise
                error('Unknown sig option "%s". Use ''none'', ''fdr'', ''bon'', or numeric alpha.', sig);
        end
    elseif isnumeric(sig) && isscalar(sig) && sig > 0 && sig < 1
        sig_method = 'fdr';   % default correction for numeric alpha
        alpha = sig;
    else
        error('sig must be ''none'', ''fdr'', ''bon'', or a numeric scalar in (0,1).');
    end

    % ---------- Basic checks ----------
    if ~ismember(varname, lme.VariableNames)
        error('Variable "%s" not found in lme.Variables.', varname);
    end

    v = lme.Variables.(varname);
    if ~iscategorical(v)
        error('Variable "%s" must be categorical.', varname);
    end

    levels   = categories(v);
    n_levels = numel(levels);
    pairs    = nchoosek(1:n_levels, 2);

    % ---------- Fixed effects: find dummy columns for varname main effect ----------
    [~, fe_tbl] = fixedEffects(lme);
    fe_names = fe_tbl.Name;                  % cellstr of fixed effect term names
    n_fixed  = numel(fe_names);

    var_prefix = [varname '_'];
    % main-effect dummy columns: start with "<varname>_" and are not interactions
    idx_cat = find(startsWith(fe_names, var_prefix) & ~contains(fe_names, ':'));

    if isempty(idx_cat)
        error('No fixed-effect dummy terms found for "%s". Check your model specification.', varname);
    end

    if numel(idx_cat) ~= (n_levels - 1)
        warning(['Expected %d dummy coefficients for "%s" but found %d. ', ...
                 'Pairwise contrasts may not map 1:1 onto category levels.'], ...
                 n_levels - 1, varname, numel(idx_cat));
    end

    % ---------- Loop over all level pairs and do Wald tests ----------
    pairwise_tbl = table();
    for i = 1:size(pairs, 1)
        i1 = pairs(i, 1);  % index in levels
        i2 = pairs(i, 2);

        c = zeros(1, n_fixed);  % contrast vector over all fixed effects

        % Default dummy coding:
        % Level 1 is reference. Coeff for level k>1 is (mean_k - mean_ref).
        %
        % Mean(level 1) = b0
        % Mean(level k>1) = b0 + b_k
        %
        % So:
        %   diff = Mean(i1) - Mean(i2)
        %   i1=1, i2>1: diff = -b_{i2}
        %   i1>1, i2=1: diff =  b_{i1}
        %   i1>1, i2>1: diff =  b_{i1} - b_{i2}

        if i1 > 1
            if (i1 - 1) <= numel(idx_cat)
                c(idx_cat(i1 - 1)) =  1;
            end
        end
        if i2 > 1
            if (i2 - 1) <= numel(idx_cat)
                c(idx_cat(i2 - 1)) = c(idx_cat(i2 - 1)) - 1;
            end
        end

        [p, F, df1, df2] = coefTest(lme, c);

        pairwise_tbl = [pairwise_tbl; ...
            table(levels(i1), levels(i2), F, df1, df2, p, ...
            'VariableNames', {'Category1', 'Category2', 'Fstat', 'DF1', 'DF2', 'pval'})]; %#ok<AGROW>
    end

    % ---------- Multiple-comparisons corrections ----------
    pairwise_tbl.pval_fdr  = mafdr(pairwise_tbl.pval, 'BHFDR', true);
    pairwise_tbl.pval_bonf = min(pairwise_tbl.pval * height(pairwise_tbl), 1);

    % For convenience: significance flag according to requested method/alpha
    switch sig_method
        case 'none'
            pairwise_tbl.sig = pairwise_tbl.pval < alpha;
            plot_field = 'pval';
        case 'fdr'
            pairwise_tbl.sig = pairwise_tbl.pval_fdr < alpha;
            plot_field = 'pval_fdr';
        case 'bon'
            pairwise_tbl.sig = pairwise_tbl.pval_bonf < alpha;
            plot_field = 'pval_bonf';
    end

    % ---------- Build symmetric matrix of p-values for plotting ----------
    pval_mat   = pivot_pairwise_matrix(pairwise_tbl, 'Category1', 'Category2', plot_field, levels);
    sig_mat    = pivot_pairwise_matrix(pairwise_tbl, 'Category1', 'Category2', 'sig',        levels);


    % Avoid log10(0)
    p_plot = pval_mat;
    p_plot(p_plot <= 0) = eps;
    logpval_mat = -log10(p_plot);

    % ---------- Plot ----------
    h = figure;
    imagesc(logpval_mat);
    axis square;
    colorbar;
    colormap(redbluecmap);  % or any diverging colormap you prefer

    xticks(1:n_levels); xticklabels(cellstr(levels));
    yticks(1:n_levels); yticklabels(cellstr(levels));

    if ~isempty(lme.ResponseName)
        respname = lme.ResponseName;
    else
        respname = 'Response';
    end

    title(sprintf('%s: -log_{10}(p_%s) for %s', ...
        respname, upper(sig_method), varname), 'Interpreter', 'none');

    % Outline significant cells
    hold on;
    [sig_y, sig_x] = find(sig_mat > 0);
    for j = 1:numel(sig_x)
        rectangle('Position', [sig_x(j)-0.5, sig_y(j)-0.5, 1, 1], ...
                  'EdgeColor', 'y', 'LineWidth', 1.5);
    end
    hold off;

    % If user doesn’t ask for the figure handle, still leave it open
    if nargout < 1
        clear h
    end

end

function M = pivot_pairwise_matrix(tbl, col1, col2, valuecol, levels)
%PIVOT_PAIRWISE_MATRIX Convert a long pairwise table into a symmetric matrix.
%
% M = pivot_pairwise_matrix(tbl, col1, col2, valuecol, levels)
%
% tbl      : table with pairwise rows (e.g., from pairwise_contrasts)
% col1     : name of first category column (e.g., 'Category1')
% col2     : name of second category column (e.g., 'Category2')
% valuecol : name of column in tbl to pivot (e.g., 'pval_fdr', 'Fstat')
% levels   : category levels in desired order (categorical, string, or cellstr)
%
% Returns:
%   M : nLevels x nLevels symmetric matrix with NaN on the diagonal.

    % ---- Normalize levels to a cellstr in desired order ----
    if iscategorical(levels)
        lvl = cellstr(categories(levels));
    elseif isstring(levels) || iscellstr(levels)
        lvl = cellstr(levels);
    else
        error('levels must be categorical, string array, or cellstr.');
    end

    n = numel(lvl);
    M = nan(n, n);

    % ---- Build a fast lookup from level label -> index ----
    idxmap = containers.Map(lvl, 1:n);

    % ---- Normalize table columns to cellstr ----
    c1_raw = tbl.(col1);
    c2_raw = tbl.(col2);

    if iscategorical(c1_raw), c1_raw = cellstr(c1_raw); end
    if iscategorical(c2_raw), c2_raw = cellstr(c2_raw); end

    if isstring(c1_raw), c1_raw = cellstr(c1_raw); end
    if isstring(c2_raw), c2_raw = cellstr(c2_raw); end

    vals = tbl.(valuecol);

    % ---- Fill matrix symmetrically ----
    for i = 1:height(tbl)
        v = vals(i);
        if isnan(v)
            continue;
        end

        key1 = c1_raw{i};
        key2 = c2_raw{i};

        if ~isKey(idxmap, key1) || ~isKey(idxmap, key2)
            continue;   % skip rows with unexpected levels
        end

        r = idxmap(key1);
        c = idxmap(key2);

        M(r, c) = v;
        M(c, r) = v;    % enforce symmetry
    end

    % Diagonal: NaN so it doesn’t clutter heatmaps
    M(1:n+1:end) = NaN;
end