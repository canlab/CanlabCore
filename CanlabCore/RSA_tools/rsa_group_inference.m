function T = rsa_group_inference(data, varargin)
% rsa_group_inference  Group-level inference on per-subject RSA scalars.
%
% Standalone Fisher-z + paired/one-sample t-test + multiple-comparisons
% correction for raw per-subject values. Useful when you already have
% per-subject scalars (not an rsm) and want the same stats machinery as
% rsm.ttest_contrasts.
%
% Two input shapes
% ----------------
%   data is [n_subjects x n_conditions] -- runs one-sample t-tests per
%                                          condition (each column independent)
%
%   data is a cell array of [n_subjects x 1] vectors -- one cell per condition
%
% Usage examples
% --------------
%   % Per-condition one-sample test against zero
%   T = rsa_group_inference(per_subject_matrix, 'names', {'hot','warm','imag'});
%
%   % With paired contrasts
%   T = rsa_group_inference(per_subject_matrix, 'names', {'HW','HI','IW'}, ...
%       'contrasts', {{'HW','HI'}; {'HW','IW'}});
%
%   % With Fisher-z transform on raw correlations
%   T = rsa_group_inference(per_subject_matrix, 'transform','fisherz');
%
% Optional name-value
% -------------------
%   'names'       cellstr of column labels (default: condition_1, ...)
%   'transform'   'none' (default) | 'fisherz'
%   'tail'        'both' (default) | 'right' | 'left'
%   'correction'  'fdr' (default) | 'bonferroni' | 'none'
%   'q'           threshold for sig flag (default 0.05)
%   'contrasts'   cell of {a,b} pairs naming named columns; each runs paired
%                 t-test of column-a minus column-b
%
% Output
% ------
%   T  table with one row per test:
%        Name, Mean, SE, t, df, P, FDR_P (or Bonf_P), sig, Cohens_d

p = inputParser;
p.addParameter('names',      {},      @(x) iscellstr(x) || isstring(x) || isempty(x)); %#ok<ISCLSTR>
p.addParameter('transform',  'none',  @(x) ischar(x) || isstring(x));
p.addParameter('tail',       'both',  @(x) ischar(x) || isstring(x));
p.addParameter('correction', 'fdr',   @(x) ischar(x) || isstring(x));
p.addParameter('q',          0.05,    @isnumeric);
p.addParameter('contrasts',  {},      @iscell);
p.parse(varargin{:});
opt = p.Results;

% Coerce data shape
if iscell(data)
    % {n_cond x 1} of [n_subj x 1] vectors
    n_cond = numel(data);
    n_subj = numel(data{1});
    X = nan(n_subj, n_cond);
    for c = 1:n_cond, X(:, c) = data{c}(:); end
elseif ismatrix(data)
    X = data;
    n_subj = size(X, 1);
    n_cond = size(X, 2);
else
    error('rsa_group_inference:badData', ...
        'data must be a matrix [n_subjects x n_conditions] or a cell of vectors.');
end

names = cellstr(opt.names);
if isempty(names)
    names = arrayfun(@(i) sprintf('condition_%d', i), 1:n_cond, 'UniformOutput', false);
elseif numel(names) ~= n_cond
    error('rsa_group_inference:badNames', ...
        'numel(names) = %d but data has %d columns.', numel(names), n_cond);
end

% Apply transform
if strcmpi(opt.transform, 'fisherz')
    X(X >  0.9999999) =  0.9999999;
    X(X < -0.9999999) = -0.9999999;
    X = atanh(X);
end

% One-sample tests per column
rows = cell(0, 1);
for c = 1:n_cond
    v = X(:, c);
    rows{end+1} = run_one_test(names{c}, v, [], opt); %#ok<AGROW>
end

% Paired contrasts
for k = 1:numel(opt.contrasts)
    pair = opt.contrasts{k};
    if numel(pair) ~= 2
        error('rsa_group_inference:badContrast', ...
            'each contrasts entry must be {a, b}.');
    end
    a_name = char(pair{1}); b_name = char(pair{2});
    ai = find(strcmp(names, a_name)); bi = find(strcmp(names, b_name));
    if isempty(ai) || isempty(bi)
        error('rsa_group_inference:contrastNotFound', ...
            'contrast {%s, %s} references unknown column.', a_name, b_name);
    end
    contrast_name = sprintf('%s_vs_%s', a_name, b_name);
    rows{end+1} = run_one_test(contrast_name, X(:, ai), X(:, bi), opt); %#ok<AGROW>
end

T = vertcat(rows{:});
T = apply_correction(T, opt);

end


function row = run_one_test(name, a, b, opt)
% Run one t-test (one-sample if b empty, paired otherwise)
if isempty(b)
    v = a;
    [~, p_i, ~, st] = ttest(v, 0, 'Tail', char(opt.tail));
    m = mean(v, 'omitnan');
    se = std(v, 0, 'omitnan') / sqrt(sum(~isnan(v)));
    d = m / std(v, 0, 'omitnan');
else
    v = a - b;
    [~, p_i, ~, st] = ttest(v, 0, 'Tail', char(opt.tail));
    m = mean(v, 'omitnan');
    se = std(v, 0, 'omitnan') / sqrt(sum(~isnan(v)));
    d = m / std(v, 0, 'omitnan');
end
row = table({name}, m, se, st.tstat, st.df, p_i, NaN, false, d, ...
    'VariableNames', {'Name', 'Mean', 'SE', 't', 'df', 'P', 'Corrected_P', 'sig', 'Cohens_d'});
end


function T = apply_correction(T, opt)
ps = T.P;
n = numel(ps);
switch lower(char(opt.correction))
    case 'fdr'
        ranks = zeros(n, 1); [~, idx] = sort(ps); ranks(idx) = 1:n;
        T.Corrected_P = min(ps(:) .* n ./ ranks, 1);
        % BH significant set
        if exist('FDR', 'file') == 2
            try, thresh = FDR(ps, opt.q); if isempty(thresh)||isnan(thresh), thresh = 0; end
                T.sig = ps <= thresh;
            catch
                T.sig = ps <= opt.q;
            end
        else
            T.sig = ps <= opt.q;
        end
        T.Properties.VariableNames{'Corrected_P'} = 'FDR_P';
    case 'bonferroni'
        T.Corrected_P = min(ps * n, 1);
        T.sig = T.Corrected_P < opt.q;
        T.Properties.VariableNames{'Corrected_P'} = 'Bonf_P';
    case 'none'
        T.Corrected_P = ps;
        T.sig = ps < opt.q;
end
end
