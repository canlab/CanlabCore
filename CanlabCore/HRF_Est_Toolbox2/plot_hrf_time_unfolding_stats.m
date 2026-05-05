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
title(sprintf('%s: %s - %s', upper(stats.model), stats.conditionA, stats.conditionB), 'Interpreter', 'none');

ax(2) = subplot(2,1,2); hold on;
plot(x, stats.p_value, 'k-');
yline(stats.alpha, 'r--', 'Alpha');
stem(x(stats.significant), stats.p_value(stats.significant), 'g.');
ylabel('p-value'); xlabel('Time bin');
title('Time-unfolding significance');
end
