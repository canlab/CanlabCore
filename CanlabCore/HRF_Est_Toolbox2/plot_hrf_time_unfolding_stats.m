function ax = plot_hrf_time_unfolding_stats(stats)
%PLOT_HRF_TIME_UNFOLDING_STATS Plot mean +/- SEM and significance mask.
figure;
ax(1) = subplot(2,1,1); hold on;
x = stats.time;
y = stats.mean;
se = stats.sem;
fill([x; flipud(x)], [y+se; flipud(y-se)], [0.6 0.6 0.9], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(x, y, 'b-', 'LineWidth', 1.8);
hline(0, 'k-');
sig_txt = '';
if isfield(stats, 'signature') && ~isempty(stats.signature)
    sig_txt = sprintf(' | %s', stats.signature);
end
unit_txt = '';
if isfield(stats, 'unit')
    unit_txt = sprintf(' | unit=%s, n=%d', stats.unit, stats.n_subjects);
end
title(sprintf('%s: %s - %s%s%s', upper(stats.model), stats.conditionA, stats.conditionB, sig_txt, unit_txt), 'Interpreter', 'none');

ax(2) = subplot(2,1,2); hold on;
plot(x, stats.p_value, 'k-');
yline(stats.alpha, 'r--', 'Alpha');
sig_idx = stats.significant & ~isnan(stats.p_value);
stem(x(sig_idx), stats.p_value(sig_idx), 'g.');
ylabel('p-value'); xlabel('Time bin');
title('Time-unfolding significance');
end
