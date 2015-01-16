function xval_lasso_brain_permutation_histogram(stats)
% Plot histograms and get permutation test-based p-value for xval_lasso_brain output structure
%
% xval_lasso_brain_permutation_histogram(stats)

create_figure('null', 2, 2);

for i = 1:4
    
    h(i) = subplot(2, 2, i);
    set(h(i), 'FontSize', 24);
    
end

% copy key info
pe = stats.full_model.pred_err;
r = stats.full_model.r_each_subject;

penull = stats.full_model.permuted.pe_values;

% permutations and p-value
nperms = length(penull);
nbins = ceil(nperms ./ 15);

pval = sum(penull <= pe) ./ nperms;
pval = max(pval, 1./nperms);

% for text on graph
xloc = mean(stats.full_model.permuted.v1000_pe_values) + std(stats.full_model.permuted.v1000_pe_values);
hi = hist(stats.full_model.permuted.v1000_pe_values, nbins);
yloc = max(hi);

% plots

subplot(2, 2, 1);
hist(stats.full_model.permuted.v1000_pe_values, nbins);
title('Perm 1000vox'); ylabel('PE');
lh = plot_vertical_line(pe, 'r');
set(lh, 'LineWidth', 3);
axis tight

text(xloc, yloc, sprintf('p < %3.6f', pval), 'FontSize', 18);

subplot(2, 2, 2);
hist(stats.full_model.permuted.v1000_r_values, nbins);
title('Perm 1000vox'); ylabel('Pred-outcome r');
lh = plot_vertical_line(r, 'r');
set(lh, 'LineWidth', 3);
axis tight

subplot(2, 2, 3);
hist(penull, nbins);
title('Perm test'); ylabel('PE');
lh = plot_vertical_line(pe, 'r');
set(lh, 'LineWidth', 3);
axis tight

subplot(2, 2, 4);
hist(stats.full_model.permuted.r_values, nbins);
title('Perm test'); ylabel('Pred-outcome r');
lh = plot_vertical_line(r, 'r');
set(lh, 'LineWidth', 3);
axis tight

equalize_axes(h([1 3]))
equalize_axes(h([2 4]))

hh = findobj(gcf, 'Type', 'patch'); set(hh, 'FaceColor', [.5 .5 .5]);

end

