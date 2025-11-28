function [h, pairwise_tbl] = pairwise_contrasts(mdl_or_stats, varname, sig)
%PAIRWISE_CONTRASTS Pairwise tests for a categorical predictor.
%
% [h, pairwise_tbl] = pairwise_contrasts(mdl_or_stats, varname, sig)
%
% mdl_or_stats : one of
%   - LinearMixedModel          (fitlme)
%   - LinearModel               (fitlm)
%   - GeneralizedLinearModel    (fitglm)
%   - ANOVA stats struct from anova1 / anovan (3rd output, often "stats")
%
% varname     : name of the categorical variable to compare (e.g., 'body_site').
%               Ignored if mdl_or_stats is an ANOVA stats struct.
%
% sig         : significance / correction option:
%                 'none'        -> uncorrected p-values, alpha = 0.05
%                 'fdr'         -> FDR-corrected p-values (Bh-FDR), alpha = 0.05 (default)
%                 'bon'         -> Bonferroni-corrected p-values, alpha = 0.05
%                 'holm_sidak'  -> Holm-Šídák FWER correction, alpha = 0.05
%                 numeric alpha -> custom alpha (FDR by default)
%
% Assumes default reference (dummy) coding for varname:
%   - First category is the reference level
%   - Fixed-effects names include varname_<Level> for levels 2..K
%
% Returns:
%   h            : figure handle for heatmap
%   pairwise_tbl : table with Category1, Category2, Mean1, Mean2, MeanDiff,
%                  Fstat, DF1, DF2, pval, pval_fdr, pval_bonf, sig,
%                  and Holm–Šídák info if requested.
%
%   Michael Sun, Ph.D. 11/26/2025 Happy Thanksgiving-Eve

    % ---------- Dispatch ANOVA stats vs model object ----------
    if isstruct(mdl_or_stats) && isfield(mdl_or_stats,'gnames')
        [h, pairwise_tbl] = pairwise_contrasts_from_anova_stats(mdl_or_stats, sig);
        return;
    end

    mdl = mdl_or_stats;

    % ---------- Model type checks ----------
    isLME = isa(mdl, 'LinearMixedModel');
    isLM  = isa(mdl, 'LinearModel');
    isGLM = isa(mdl, 'GeneralizedLinearModel');

    if ~(isLME || isLM || isGLM)
        error(['mdl_or_stats must be a LinearModel, LinearMixedModel, ', ...
               'GeneralizedLinearModel, or ANOVA stats struct.']);
    end

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
                sig_method = 'none';        alpha = 0.05;
            case "fdr"
                sig_method = 'fdr';         alpha = 0.05;
            case {"bon","bonf","bonferroni"}
                sig_method = 'bon';         alpha = 0.05;
            case "holm_sidak"
                sig_method = 'holm_sidak';  alpha = 0.05;
            otherwise
                error('Unknown sig option "%s". Use ''none'', ''fdr'', ''bon'', ''holm_sidak'', or numeric alpha.',char(sig));
        end
    elseif isnumeric(sig) && isscalar(sig) && sig > 0 && sig < 1
        sig_method = 'fdr';
        alpha = sig;
    else
        error('sig must be ''none'', ''fdr'', ''bon'', ''holm_sidak'', or a numeric scalar in (0,1).');
    end

    % ---------- Basic checks ----------
    if ~ismember(varname, mdl.VariableNames)
        error('Variable "%s" not found in mdl.Variables.', varname);
    end

    v = mdl.Variables.(varname);
    if ~iscategorical(v)
        error('Variable "%s" must be categorical.', varname);
    end

    levels   = categories(v);
    n_levels = numel(levels);
    pairs    = nchoosek(1:n_levels, 2);

    % ---------- Fixed-effect coefficients & names ----------
    if isLME
        [beta, fe_tbl] = fixedEffects(mdl);
        fe_names = fe_tbl.Name;
    else
        beta    = mdl.Coefficients.Estimate;
        fe_names = mdl.CoefficientNames(:);
    end
    n_fixed = numel(fe_names);

    % ---------- Find intercept index ----------
    idx_int = find(strcmp(fe_names,'(Intercept)'));
    if isempty(idx_int)
        warning('No "(Intercept)" term found in fixed effects. Using baseline = 0.');
        baseline = 0;
    else
        baseline = beta(idx_int);
    end

    % ---------- Find dummy coeffs for varname (reference coding assumed) ----------
    % Names like 'body_site_Right Face', 'body_site_Left Arm', ...
    var_prefix = [varname '_'];
    idx_cat = find(startsWith(fe_names, var_prefix) & ~contains(fe_names,':'));

    if isempty(idx_cat)
        error('No fixed-effect dummy terms found for "%s". Check your model specification/coding.', varname);
    end

    if numel(idx_cat) ~= (n_levels - 1)
        warning(['Expected %d dummy coefficients for "%s" but found %d. ', ...
                 'The mapping may be off if you are not using default reference coding.'], ...
                 n_levels - 1, varname, numel(idx_cat));
    end

    % For sanity, print mapping if you like:
    % disp(table(levels(2:end), fe_names(idx_cat),'VariableNames',{'Level','CoeffName'}));

    % ---------- Compute level means on the linear scale ----------
    % Reference level (level 1) is baseline.
    level_means = NaN(n_levels,1);
    level_means(1) = baseline;

    for k = 2:n_levels
        if (k-1) <= numel(idx_cat)
            level_means(k) = baseline + beta(idx_cat(k-1));
        else
            level_means(k) = NaN;
        end
    end

    % ---------- Loop over all level pairs and build contrasts ----------
    pairwise_tbl = table();
    tol = 1e-12;

    for i = 1:size(pairs,1)
        i1 = pairs(i,1);
        i2 = pairs(i,2);

        c = zeros(1, n_fixed);  % contrast over fixed effects

        % Reference/dummy logic:
        %   Level 1: mean = baseline
        %   Level k>1: mean = baseline + beta_cat(k-1)
        %
        % diff = mean(i1) - mean(i2)
        %   i1=1, i2>1: diff = -beta_cat(i2-1)
        %   i1>1, i2=1: diff =  beta_cat(i1-1)
        %   i1>1, i2>1: diff =  beta_cat(i1-1) - beta_cat(i2-1)

        if i1 > 1 && (i1-1) <= numel(idx_cat)
            c(idx_cat(i1-1)) = 1;
        end
        if i2 > 1 && (i2-1) <= numel(idx_cat)
            c(idx_cat(i2-1)) = c(idx_cat(i2-1)) - 1;
        end

        Mean1    = level_means(i1);
        Mean2    = level_means(i2);
        MeanDiff = Mean1 - Mean2;   % Category1 - Category2

        if all(abs(c) < tol)
            % No testable contrast (should not happen with proper coding, but be safe)
            p   = NaN; F = NaN; df1 = NaN; df2 = NaN;
        else
            [p, F, df1, df2] = coefTest(mdl, c);
        end

        pairwise_tbl = [pairwise_tbl; ...
            table(levels(i1), levels(i2), ...
                  Mean1, Mean2, MeanDiff, ...
                  F, df1, df2, p, ...
            'VariableNames', {'Category1','Category2', ...
                              'Mean1','Mean2','MeanDiff', ...
                              'Fstat','DF1','DF2','pval'})]; %#ok<AGROW>
    end

    % ---------- Corrections + plotting ----------
    [h, pairwise_tbl] = add_corrections_and_plot(mdl, pairwise_tbl, levels, varname, sig_method, alpha);
end

function [h, pairwise_tbl] = pairwise_contrasts_from_anova_stats(stats, sig)

    % ----- sig handling -----
    if nargin < 2 || isempty(sig)
        sig_method = 'fdr'; alpha = 0.05;
    elseif ischar(sig) || isstring(sig)
        sig = lower(string(sig));
        switch sig
            case "none",       sig_method = 'none';       alpha = 0.05;
            case "fdr",        sig_method = 'fdr';        alpha = 0.05;
            case {"bon","bonf","bonferroni"}
                              sig_method = 'bon';        alpha = 0.05;
            case "holm_sidak", sig_method = 'holm_sidak'; alpha = 0.05;
            otherwise
                error('Unknown sig option "%s".', char(sig));
        end
    elseif isnumeric(sig) && isscalar(sig) && sig > 0 && sig < 1
        sig_method = 'fdr'; alpha = sig;
    else
        error('sig must be ''none'', ''fdr'', ''bon'', ''holm_sidak'', or a numeric scalar.');
    end

    levels = stats.gnames(:);

    % multcompare: c = [i1,i2,lowerCI,diff,upperCI,pval]
    [c, ~, ~, gnames] = multcompare(stats, 'alpha', alpha, 'display','off');

    Category1 = gnames(c(:,1));
    Category2 = gnames(c(:,2));
    pvals     = c(:,6);

    if isfield(stats,'means')
        gm = stats.means(:);
        Mean1    = gm(c(:,1));
        Mean2    = gm(c(:,2));
        MeanDiff = Mean1 - Mean2;
    else
        Mean1 = NaN(size(pvals));
        Mean2 = NaN(size(pvals));
        MeanDiff = NaN(size(pvals));
    end

    if isfield(stats,'df') && numel(stats.df) >= 2
        df_error = stats.df(2);
    else
        df_error = NaN;
    end

    Fstat = NaN(size(pvals));
    DF1   = ones(size(pvals));
    DF2   = repmat(df_error, size(pvals));

    pairwise_tbl = table(Category1, Category2, ...
                         Mean1, Mean2, MeanDiff, ...
                         Fstat, DF1, DF2, pvals, ...
        'VariableNames', {'Category1','Category2', ...
                          'Mean1','Mean2','MeanDiff', ...
                          'Fstat','DF1','DF2','pval'});

    mdl_placeholder.ResponseName  = 'Response';
    mdl_placeholder.VariableNames = {'Factor'};

    [h, pairwise_tbl] = add_corrections_and_plot(mdl_placeholder, pairwise_tbl, levels, 'Factor', sig_method, alpha);
end

function [h, pairwise_tbl] = add_corrections_and_plot(mdl, pairwise_tbl, levels, varname, sig_method, alpha)

    % ----- Multiple-comparisons corrections -----
    pairwise_tbl.pval_fdr  = mafdr(pairwise_tbl.pval, 'BHFDR', true);
    pairwise_tbl.pval_bonf = min(pairwise_tbl.pval * height(pairwise_tbl), 1);

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

        case 'holm_sidak'
            if exist('holm_sidak','file') ~= 2
                error('sig option ''holm_sidak'' requires holm_sidak.m on the path.');
            end
            [sig_hs, pthr_hs] = holm_sidak(pairwise_tbl.pval, alpha);
            pairwise_tbl.sig_holm_sidak  = sig_hs(:);
            pairwise_tbl.holm_sidak_pthr = repmat(pthr_hs, height(pairwise_tbl), 1);
            pairwise_tbl.sig = pairwise_tbl.sig_holm_sidak;
            plot_field       = 'pval';  % still plotting raw p
    end

    % ----- Build matrices -----
    pval_mat = pivot_pairwise_matrix(pairwise_tbl, 'Category1', 'Category2', plot_field, levels);
    sig_mat  = pivot_pairwise_matrix(pairwise_tbl, 'Category1', 'Category2', 'sig',       levels);

    p_plot = pval_mat;
    p_plot(p_plot <= 0) = eps;
    logpval_mat = -log10(p_plot);

    % ----- Plot -----
    if iscategorical(levels)
        lvl_labels = cellstr(levels);
    elseif isstring(levels) || iscellstr(levels)
        lvl_labels = cellstr(levels);
    else
        error('levels must be categorical, string array, or cellstr.');
    end

    n = numel(lvl_labels);

    h = figure;
    imagesc(logpval_mat);
    axis square;
    colorbar;
    colormap(redbluecmap);

    xticks(1:n); xticklabels(lvl_labels);
    yticks(1:n); yticklabels(lvl_labels);

    if isfield(mdl,'ResponseName') && ~isempty(mdl.ResponseName)
        respname = mdl.ResponseName;
    else
        respname = 'Response';
    end

    title(sprintf('%s: -log_{10}(p_%s) for %s', ...
        respname, upper(sig_method), varname), 'Interpreter','none');

    hold on;
    [sy, sx] = find(sig_mat > 0);
    for j = 1:numel(sx)
        rectangle('Position',[sx(j)-0.5, sy(j)-0.5, 1, 1], ...
                  'EdgeColor','y','LineWidth',1.5);
    end
    hold off;
end

function M = pivot_pairwise_matrix(tbl, col1, col2, valuecol, levels)

    if iscategorical(levels)
        lvl = cellstr(levels);
    elseif isstring(levels) || iscellstr(levels)
        lvl = cellstr(levels);
    else
        error('levels must be categorical, string array, or cellstr.');
    end

    n = numel(lvl);
    M = nan(n, n);

    idxmap = containers.Map(lvl, 1:n);

    c1_raw = tbl.(col1);
    c2_raw = tbl.(col2);

    if iscategorical(c1_raw), c1_raw = cellstr(c1_raw); end
    if iscategorical(c2_raw), c2_raw = cellstr(c2_raw); end
    if isstring(c1_raw),     c1_raw = cellstr(c1_raw); end
    if isstring(c2_raw),     c2_raw = cellstr(c2_raw); end

    vals = tbl.(valuecol);

    for i = 1:height(tbl)
        v = vals(i);
        if isnan(v), continue; end

        key1 = c1_raw{i};
        key2 = c2_raw{i};

        if ~isKey(idxmap,key1) || ~isKey(idxmap,key2), continue; end

        r = idxmap(key1);
        c = idxmap(key2);

        M(r, c) = v;
        M(c, r) = v;
    end

    M(1:n+1:end) = NaN;
end