function stats = hrf_2x2_study_score_stats(a1b1, a1b2, a2b1, a2b2, varargin)
%HRF_2X2_STUDY_SCORE_STATS Repeated-measures 2x2 tests for study scores.
%
% stats = hrf_2x2_study_score_stats(a1b1, a1b2, a2b1, a2b2, ...)
%
% The four study inputs should be study structs returned by
% hrf_input_table_to_study or hrf_second_level_inputs_to_study. Each cell is
% reduced to one curve per subject using hrf_time_unfolding_stats, subjects
% are aligned across all four cells, and paired tests are run for simple
% effects, main effects, and the 2x2 interaction.
%
% Cell layout:
%   a1b1  a1b2
%   a2b1  a2b2
%
% Example:
% stats = hrf_2x2_study_score_stats(study_scores_lf_exp, study_scores_lf_acc, ...
%     study_scores_obs_exp, study_scores_obs_acc, ...
%     'FactorA', {'lf','obs'}, 'FactorB', {'exp','acc'}, ...
%     'Model', 'sfir', 'Signature', 'sig_all_NPS', ...
%     'Condition', '*heat_start_ttl_1');

p = inputParser;
p.addRequired('a1b1', @isstruct);
p.addRequired('a1b2', @isstruct);
p.addRequired('a2b1', @isstruct);
p.addRequired('a2b2', @isstruct);
p.addParameter('FactorA', {'A1', 'A2'}, @(x) iscell(x) || isstring(x));
p.addParameter('FactorB', {'B1', 'B2'}, @(x) iscell(x) || isstring(x));
p.addParameter('FactorAName', 'Factor A', @(x) ischar(x) || isstring(x));
p.addParameter('FactorBName', 'Factor B', @(x) ischar(x) || isstring(x));
p.addParameter('CellLabels', {}, @(x) iscell(x) || isstring(x));
p.addParameter('Model', 'sfir', @(x) ischar(x) || isstring(x));
p.addParameter('SourceModel', '', @(x) ischar(x) || isstring(x));
p.addParameter('Signature', '', @(x) ischar(x) || isstring(x));
p.addParameter('Condition', '', @(x) ischar(x) || isstring(x));
p.addParameter('ConditionA', '', @(x) ischar(x) || isstring(x));
p.addParameter('ConditionB', '', @(x) ischar(x) || isstring(x));
p.addParameter('Unit', 'subject', @(x) ischar(x) || isstring(x));
p.addParameter('MissingPolicy', 'warn', @(x) ischar(x) || isstring(x));
p.addParameter('Alpha', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
% Across-lag correction for the 2x2 contrast timecourses (shared engine):
% 'none' (default) | 'fdr' | 'permutation' | 'cluster'. When not 'none', each
% contrast gains .p_corrected / .significant_corrected.
p.addParameter('Correction', 'none', @(x) ischar(x) || isstring(x));
p.addParameter('Nperm', 5000, @(x) isscalar(x) && x >= 100);
p.parse(a1b1, a1b2, a2b1, a2b2, varargin{:});
opts = p.Results;

factorA = local_two_names(opts.FactorA, 'FactorA');
factorB = local_two_names(opts.FactorB, 'FactorB');
conditionA = char(opts.ConditionA);
if isempty(conditionA)
    conditionA = char(opts.Condition);
end
unit = strtrim(char(opts.Unit));
if ~strcmpi(unit, 'subject')
    error('hrf_2x2_study_score_stats currently supports Unit=''subject'' so paired tests align one curve per subject.');
end

cell_labels = local_cell_labels(opts.CellLabels, factorA, factorB);
source_stats = cell(2, 2);
source_stats{1, 1} = local_cell_time_stats(a1b1, opts, conditionA);
source_stats{1, 2} = local_cell_time_stats(a1b2, opts, conditionA);
source_stats{2, 1} = local_cell_time_stats(a2b1, opts, conditionA);
source_stats{2, 2} = local_cell_time_stats(a2b2, opts, conditionA);

subject_ids = local_common_subjects(source_stats);
time = local_common_time(source_stats);

Y11 = local_align_data(source_stats{1, 1}, subject_ids);
Y12 = local_align_data(source_stats{1, 2}, subject_ids);
Y21 = local_align_data(source_stats{2, 1}, subject_ids);
Y22 = local_align_data(source_stats{2, 2}, subject_ids);

stats = struct();
stats.time = time(:);
stats.alpha = opts.Alpha;
stats.model = char(opts.Model);
stats.source_model = char(opts.SourceModel);
stats.signature = char(opts.Signature);
stats.conditionA = conditionA;
stats.conditionB = char(opts.ConditionB);
stats.unit = char(opts.Unit);
stats.factorA = factorA(:);
stats.factorB = factorB(:);
stats.factorA_name = char(opts.FactorAName);
stats.factorB_name = char(opts.FactorBName);
stats.subject_ids = subject_ids(:);
stats.n_subjects = numel(subject_ids);

stats.cells = struct();
stats.cells.A1B1 = local_cell_summary(cell_labels{1}, factorA{1}, factorB{1}, Y11, source_stats{1, 1});
stats.cells.A1B2 = local_cell_summary(cell_labels{2}, factorA{1}, factorB{2}, Y12, source_stats{1, 2});
stats.cells.A2B1 = local_cell_summary(cell_labels{3}, factorA{2}, factorB{1}, Y21, source_stats{2, 1});
stats.cells.A2B2 = local_cell_summary(cell_labels{4}, factorA{2}, factorB{2}, Y22, source_stats{2, 2});

contrasts = struct();
contrasts.A1_B1_minus_B2 = local_contrast_stats( ...
    sprintf('%s: %s - %s', factorA{1}, factorB{1}, factorB{2}), ...
    'A1_B1_minus_B2', 'A1B1 - A1B2', Y11 - Y12, opts.Alpha);
contrasts.A2_B1_minus_B2 = local_contrast_stats( ...
    sprintf('%s: %s - %s', factorA{2}, factorB{1}, factorB{2}), ...
    'A2_B1_minus_B2', 'A2B1 - A2B2', Y21 - Y22, opts.Alpha);
contrasts.B1_A1_minus_A2 = local_contrast_stats( ...
    sprintf('%s: %s - %s', factorB{1}, factorA{1}, factorA{2}), ...
    'B1_A1_minus_A2', 'A1B1 - A2B1', Y11 - Y21, opts.Alpha);
contrasts.B2_A1_minus_A2 = local_contrast_stats( ...
    sprintf('%s: %s - %s', factorB{2}, factorA{1}, factorA{2}), ...
    'B2_A1_minus_A2', 'A1B2 - A2B2', Y12 - Y22, opts.Alpha);
contrasts.main_A = local_contrast_stats( ...
    sprintf('Main %s: %s - %s', stats.factorA_name, factorA{1}, factorA{2}), ...
    'main_A', 'mean(A1B1,A1B2) - mean(A2B1,A2B2)', ...
    (Y11 + Y12) ./ 2 - (Y21 + Y22) ./ 2, opts.Alpha);
contrasts.main_B = local_contrast_stats( ...
    sprintf('Main %s: %s - %s', stats.factorB_name, factorB{1}, factorB{2}), ...
    'main_B', 'mean(A1B1,A2B1) - mean(A1B2,A2B2)', ...
    (Y11 + Y21) ./ 2 - (Y12 + Y22) ./ 2, opts.Alpha);
contrasts.interaction_AxB = local_contrast_stats( ...
    sprintf('Interaction: (%s - %s) x (%s - %s)', factorA{1}, factorA{2}, factorB{1}, factorB{2}), ...
    'interaction_AxB', '(A1B1 - A1B2) - (A2B1 - A2B2)', ...
    (Y11 - Y12) - (Y21 - Y22), opts.Alpha);

% Across-lag correction of every contrast timecourse via the shared engine
% (each contrast is a one-sample test of its paired difference vs 0). Additive:
% the uncorrected per-timepoint .significant is preserved.
stats.correction = lower(char(opts.Correction));
if ~strcmpi(stats.correction, 'none')
    fn = fieldnames(contrasts);
    for i = 1:numel(fn)
        Cc = hrf_time_correction(contrasts.(fn{i}).data, 'Correction', stats.correction, ...
            'Nperm', opts.Nperm, 'Alpha', opts.Alpha);
        contrasts.(fn{i}).p_corrected = Cc.p_corr(:);
        contrasts.(fn{i}).significant_corrected = Cc.sig(:);
        contrasts.(fn{i}).correction = stats.correction;
    end
end

stats.contrasts = contrasts;
stats.contrast_names = fieldnames(contrasts);
stats.contrast_table = local_contrast_table(contrasts);
stats.source_stats = source_stats;
end

function cell_stats = local_cell_time_stats(study, opts, conditionA)
cell_stats = hrf_time_unfolding_stats(study, ...
    'Model', opts.Model, ...
    'SourceModel', opts.SourceModel, ...
    'Signature', opts.Signature, ...
    'ConditionA', conditionA, ...
    'ConditionB', opts.ConditionB, ...
    'Unit', opts.Unit, ...
    'MissingPolicy', opts.MissingPolicy, ...
    'Alpha', opts.Alpha);
end

function names = local_two_names(input_names, param_name)
names = cellstr(string(input_names));
names = names(:)';
if numel(names) ~= 2 || any(cellfun(@isempty, names))
    error('%s must contain exactly two non-empty names.', param_name);
end
end

function labels = local_cell_labels(input_labels, factorA, factorB)
if isempty(input_labels)
    labels = { ...
        sprintf('%s_%s', factorA{1}, factorB{1}), ...
        sprintf('%s_%s', factorA{1}, factorB{2}), ...
        sprintf('%s_%s', factorA{2}, factorB{1}), ...
        sprintf('%s_%s', factorA{2}, factorB{2})};
    return
end
labels = cellstr(string(input_labels));
labels = labels(:)';
if numel(labels) ~= 4 || any(cellfun(@isempty, labels))
    error('CellLabels must contain exactly four non-empty labels.');
end
end

function subject_ids = local_common_subjects(source_stats)
subject_ids = cellstr(string(source_stats{1}.subject_ids(:)));
for i = 2:numel(source_stats)
    this_subjects = cellstr(string(source_stats{i}.subject_ids(:)));
    keep = ismember(subject_ids, this_subjects);
    subject_ids = subject_ids(keep);
end
if isempty(subject_ids)
    error('No common subjects were found across the four study cells.');
end
end

function time = local_common_time(source_stats)
time = source_stats{1}.time(:);
for i = 2:numel(source_stats)
    this_time = source_stats{i}.time(:);
    if numel(this_time) ~= numel(time) || any(abs(this_time - time) > max(eps(time)))
        error('Time axes do not match across the four study cells.');
    end
end
end

function Y = local_align_data(cell_stats, subject_ids)
cell_subjects = cellstr(string(cell_stats.subject_ids(:)));
[tf, loc] = ismember(subject_ids, cell_subjects);
if any(~tf)
    error('Subject alignment failed unexpectedly.');
end
Y = cell_stats.data(loc, :);
end

function summary = local_cell_summary(label, factorA, factorB, data, source_stats)
summary = struct();
summary.label = label;
summary.factorA = factorA;
summary.factorB = factorB;
summary.data = data;
summary.mean = local_mean_omitnan(data, 1)';
summary.sem = local_sem_omitnan(data, 1)';
summary.n_subjects = sum(~isnan(data), 1)';
summary.source_stats = source_stats;
end

function c = local_contrast_stats(name, fieldname, formula, data, alpha)
[T, P] = local_ttest_each_timepoint(data);
c = struct();
c.name = name;
c.fieldname = fieldname;
c.formula = formula;
c.data = data;
c.mean = local_mean_omitnan(data, 1)';
c.sem = local_sem_omitnan(data, 1)';
c.t_value = T(:);
c.p_value = P(:);
c.significant = P(:) < alpha;
c.n_subjects = sum(~isnan(data), 1)';
c.alpha = alpha;
end

function tbl = local_contrast_table(contrasts)
names = fieldnames(contrasts);
label = cell(numel(names), 1);
formula = cell(numel(names), 1);
n_significant = zeros(numel(names), 1);
n_significant_corrected = nan(numel(names), 1);
min_p = nan(numel(names), 1);
for i = 1:numel(names)
    c = contrasts.(names{i});
    label{i} = c.name;
    formula{i} = c.formula;
    n_significant(i) = sum(c.significant & ~isnan(c.p_value));
    if isfield(c, 'significant_corrected')
        n_significant_corrected(i) = sum(c.significant_corrected);
    end
    valid_p = c.p_value(~isnan(c.p_value));
    if isempty(valid_p)
        min_p(i) = NaN;
    else
        min_p(i) = min(valid_p);
    end
end
tbl = table(names(:), label, formula, n_significant, n_significant_corrected, min_p, ...
    'VariableNames', {'contrast', 'label', 'formula', 'n_significant_timepoints', ...
    'n_significant_corrected', 'min_p'});
end

function [T, P] = local_ttest_each_timepoint(D)
n_tp = size(D, 2);
P = nan(1, n_tp);
T = nan(1, n_tp);
for t = 1:n_tp
    y = D(:, t);
    y = y(~isnan(y));
    if numel(y) < 2
        continue
    end
    [~, pval, ~, st] = ttest(y);
    P(t) = pval;
    T(t) = st.tstat;
end
end

function m = local_mean_omitnan(X, dim)
if nargin < 2
    dim = 1;
end
valid = ~isnan(X);
den = sum(valid, dim);
X(~valid) = 0;
m = sum(X, dim) ./ den;
m(den == 0) = NaN;
end

function se = local_sem_omitnan(X, dim)
if nargin < 2
    dim = 1;
end
mu = local_mean_omitnan(X, dim);
if dim == 1
    centered = X - repmat(mu, size(X, 1), 1);
elseif dim == 2
    centered = X - repmat(mu, 1, size(X, 2));
else
    error('local_sem_omitnan supports dim 1 or 2.');
end
centered(isnan(centered)) = 0;
n = sum(~isnan(X), dim);
sd = sqrt(sum(centered .^ 2, dim) ./ max(n - 1, 1));
se = sd ./ sqrt(n);
se(n == 0) = NaN;
end
