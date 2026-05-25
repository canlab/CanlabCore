function ax = plot_hrf_2x2_study_score_stats(stats, varargin)
%PLOT_HRF_2X2_STUDY_SCORE_STATS Plot one contrast from 2x2 study-score stats.

p = inputParser;
p.addRequired('stats', @isstruct);
p.addParameter('Contrast', 'interaction_AxB', @(x) ischar(x) || isstring(x));
p.addParameter('ShowCells', false, @(x) islogical(x) || isnumeric(x));
p.parse(stats, varargin{:});
opts = p.Results;

contrast_field = local_contrast_field(stats, char(opts.Contrast));
c = stats.contrasts.(contrast_field);
x = stats.time(:);

figure;
ax(1) = subplot(2, 1, 1);
hold on;
if logical(opts.ShowCells)
    local_plot_cell_means(stats, x);
end
fill([x; flipud(x)], [c.mean + c.sem; flipud(c.mean - c.sem)], ...
    [0.3 0.3 0.3], 'FaceAlpha', 0.18, 'EdgeColor', 'none');
plot(x, c.mean, 'k-', 'LineWidth', 2);
line([min(x) max(x)], [0 0], 'Color', [0 0 0], 'LineStyle', '-');
ylabel('Difference in pattern score');
title(local_title(stats, c), 'Interpreter', 'none');

ax(2) = subplot(2, 1, 2);
hold on;
plot(x, c.p_value, 'k-', 'LineWidth', 1.2);
line([min(x) max(x)], [stats.alpha stats.alpha], 'Color', [0.8 0 0], 'LineStyle', '--');
sig_idx = c.significant & ~isnan(c.p_value);
stem(x(sig_idx), c.p_value(sig_idx), 'g.');
ylabel('p-value');
xlabel('Time bin');
title(sprintf('%s, n=%d', c.formula, stats.n_subjects), 'Interpreter', 'none');
end

function field = local_contrast_field(stats, query)
fields = fieldnames(stats.contrasts);
if ismember(query, fields)
    field = query;
    return
end
for i = 1:numel(fields)
    c = stats.contrasts.(fields{i});
    if strcmpi(query, c.name) || strcmpi(query, c.formula)
        field = fields{i};
        return
    end
end
error('Unknown contrast "%s". Available contrasts: %s', query, strjoin(fields, ', '));
end

function local_plot_cell_means(stats, x)
cell_fields = fieldnames(stats.cells);
colors = lines(numel(cell_fields));
for i = 1:numel(cell_fields)
    cell_stats = stats.cells.(cell_fields{i});
    plot(x, cell_stats.mean, '-', 'Color', colors(i, :), 'LineWidth', 0.9, ...
        'DisplayName', cell_stats.label);
end
legend('Location', 'best', 'Interpreter', 'none');
end

function txt = local_title(stats, c)
sig_txt = '';
if isfield(stats, 'signature') && ~isempty(stats.signature)
    sig_txt = sprintf(' | %s', stats.signature);
end
condition_txt = stats.conditionA;
if isempty(condition_txt)
    condition_txt = 'first condition';
end
txt = sprintf('%s | %s%s', c.name, condition_txt, sig_txt);
end
