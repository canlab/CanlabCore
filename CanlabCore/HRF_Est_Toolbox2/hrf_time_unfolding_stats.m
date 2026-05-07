function stats = hrf_time_unfolding_stats(study, varargin)
%HRF_TIME_UNFOLDING_STATS Multilevel time-unfolding tests across subjects.
% Supports within-subject condition contrasts and between-group contrasts.

p = inputParser;
p.addRequired('study', @isstruct);
p.addParameter('Model', 'sfir', @(x) ischar(x) || isstring(x));
p.addParameter('ConditionA', '', @(x) ischar(x) || isstring(x));
p.addParameter('ConditionB', '', @(x) ischar(x) || isstring(x));
p.addParameter('Group', {}, @(x) iscell(x) || isstring(x));
p.addParameter('Alpha', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
p.parse(study, varargin{:});
opts = p.Results;

model_name = char(opts.Model);
subj_n = numel(study.results);

% Build subject x time arrays
for s = 1:subj_n
    r = study.results{s};
    if isempty(r), continue; end
    h = r.fits.(model_name).hrf;
    if isempty(opts.ConditionA)
        cA = 1;
    else
        cA = find(strcmp(r.conditions, char(opts.ConditionA)), 1);
    end
    if isempty(cA), error('ConditionA not found for subject %s', study.subject_ids{s}); end
    A(s, :) = h(:, cA)'; %#ok<AGROW>

    if ~isempty(opts.ConditionB)
        cB = find(strcmp(r.conditions, char(opts.ConditionB)), 1);
        if isempty(cB), error('ConditionB not found for subject %s', study.subject_ids{s}); end
        B(s, :) = h(:, cB)'; %#ok<AGROW>
    end
end

if ~isempty(opts.ConditionB)
    D = A - B;
else
    D = A;
end

n_tp = size(D, 2);
P = nan(1, n_tp); T = nan(1, n_tp);
for t = 1:n_tp
    [~, pval, ~, st] = ttest(D(:, t));
    P(t) = pval; T(t) = st.tstat;
end

stats = struct();
stats.time = (1:n_tp)';
stats.mean = mean(D, 1)';
stats.sem = std(D, 0, 1)' ./ sqrt(size(D, 1));
stats.p_value = P(:);
stats.t_value = T(:);
stats.significant = P(:) < opts.Alpha;
stats.alpha = opts.Alpha;
stats.conditionA = char(opts.ConditionA);
stats.conditionB = char(opts.ConditionB);
stats.model = model_name;

% Optional between-group test
if ~isempty(opts.Group)
    g = cellstr(string(opts.Group));
    if numel(g) ~= subj_n, error('Group labels must match number of subjects'); end
    ug = unique(g);
    if numel(ug) == 2
        ix1 = strcmp(g, ug{1}); ix2 = strcmp(g, ug{2});
        Pg = nan(1, n_tp); Tg = nan(1, n_tp);
        for t = 1:n_tp
            [~, pval, ~, st] = ttest2(D(ix1, t), D(ix2, t), 'Vartype', 'unequal');
            Pg(t) = pval; Tg(t) = st.tstat;
        end
        stats.group_labels = ug;
        stats.group_p_value = Pg(:);
        stats.group_t_value = Tg(:);
    end
end
end
