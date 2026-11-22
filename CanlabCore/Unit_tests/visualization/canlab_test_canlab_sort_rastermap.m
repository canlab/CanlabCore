function tests = canlab_test_canlab_sort_rastermap
%CANLAB_TEST_CANLAB_SORT_RASTERMAP Rastermap row sorting (canlab_sort_rastermap.m) and the
% fmri_data.plot 'rastermap' carpet option.
tests = functiontests(localfunctions);
end


function setup(~)
close all
end


function teardown(~)
close all
end


function [X, phase] = ring_data(n, t)
% Rows tuned to a phase on a ring; the stimulus drifts around the ring.
rng(5);
phase = 2*pi*rand(n, 1);
stim = cumsum(randn(1, t)) / 8;
X = exp(2*cos(phase - stim)) + randn(n, t);
end


function test_recovers_ring_order(tc)
[X, phase] = ring_data(2000, 500);
isort = canlab_sort_rastermap(X, 'noverbose');
tc.verifySize(isort, [2000 1]);
tc.verifyEqual(sort(isort), (1:2000)');
d = abs(angle(exp(1i * diff(phase(isort))))) * 180 / pi;
% random order gives ~90 deg between neighbours; the sort should be smooth
tc.verifyLessThan(mean(d), 10);
end


function test_outputs_and_invalid_rows(tc)
[X, ~] = ring_data(600, 200);
X(5, :) = NaN;                 % non-finite row
X(9, :) = 3;                   % constant row
[isort, labels, OUT] = canlab_sort_rastermap(X, 'noverbose', 'n_clusters', 30);
tc.verifyFalse(OUT.igood(5));
tc.verifyFalse(OUT.igood(9));
tc.verifyTrue(all(ismember([5 9], isort(end-1:end))));     % invalid rows last
tc.verifyTrue(all(isnan(OUT.embedding([5 9]))));
tc.verifyEmpty(labels);
tc.verifyEqual(size(OUT.X_embedding, 2), 200);
tc.verifyTrue(all(isfield(OUT.timing, {'pca', 'kmeans', 'sort', 'upsample', 'total'})));
end


function test_hierarchical_cluster_labels(tc)
[X, ~] = ring_data(1500, 300);
[isort, labels, OUT] = canlab_sort_rastermap(X, 'noverbose', 'k', 5);
tc.verifySize(labels, [1500 1]);
tc.verifyEqual(unique(labels)', 1:5);
% labels are numbered in order of first appearance down the sorted rows
ls = labels(isort);
[~, first] = unique(ls, 'stable');
tc.verifyEqual(ls(sort(first))', 1:5);
tc.verifyNotEmpty(OUT.linkage_tree);
tc.verifyNotEmpty(OUT.bin_labels);
end


function test_contiguous_clusters_are_bands(tc)
[X, ~] = ring_data(1500, 300);
[isort, labels] = canlab_sort_rastermap(X, 'noverbose', 'cluster', 4, 'cluster_method', 'contiguous');
ls = labels(isort);
% each cluster is one contiguous band: exactly k-1 label changes
tc.verifyEqual(sum(diff(ls) ~= 0), 3);
end


function test_sequence_preset_and_overrides(tc)
[X, ~] = ring_data(800, 300);
[~, ~, OUT] = canlab_sort_rastermap(X, 'noverbose', 'sequence');
tc.verifyEqual(OUT.settings.locality, 0.75);
tc.verifyEqual(OUT.settings.time_lag_window, 10);
[~, ~, OUT] = canlab_sort_rastermap(X, 'noverbose', 'sequence', 'locality', 0.5);
tc.verifyEqual(OUT.settings.locality, 0.5);
end


function test_reproducible_and_rng_restored(tc)
[X, ~] = ring_data(800, 200);
rng(123); r0 = rand;
rng(123);
i1 = canlab_sort_rastermap(X, 'noverbose');
r1 = rand;                       % caller's RNG stream continues unaffected
i2 = canlab_sort_rastermap(X, 'noverbose');
tc.verifyEqual(i1, i2);
tc.verifyEqual(r1, r0);
end


function test_plot_with_clusters_and_colors(tc)
[X, ~] = ring_data(800, 200);
[~, ~, OUT] = canlab_sort_rastermap(X, 'noverbose', 'k', 3, 'plot', 'colors', {[1 0 0], [0 1 0], [0 0 1]});
tc.verifyTrue(isgraphics(OUT.fig_handle));
[~, ~, OUT] = canlab_sort_rastermap(X, 'noverbose', 'k', 3, 'plot', 'colors', 'lines');
tc.verifyTrue(isgraphics(OUT.fig_handle));
[~, ~, OUT] = canlab_sort_rastermap(X, 'noverbose', 'plot');
tc.verifyTrue(isgraphics(OUT.fig_handle));
end


function test_few_rows_skip_clustering(tc)
rng(2);
X = randn(40, 60);
tc.verifyWarning(@() canlab_sort_rastermap(X, 'noverbose'), 'canlab_sort_rastermap:fewRows');
end


function test_fmri_data_plot_rastermap(tc)
obj = canlab_get_sample_fmri_data();
[~, ~, Xout] = plot(obj, 'rastermap', 'noorthviews', 'nomontages');
tc.verifyTrue(isfield(Xout, 'voxel_sort_order'));
tc.verifyEqual(sort(Xout.voxel_sort_order(:)), (1:size(obj.dat, 1))');
end


function test_plot_correlation_matrix_reorder_by_rastermap(tc)
% 3 blocks of correlated variables, shuffled: the reorder should group them
rng(4);
nobs = 200; blocks = repelem(1:3, 8); latent = randn(nobs, 3);
X = latent(:, blocks) + 0.7 * randn(nobs, numel(blocks));
shuf = randperm(numel(blocks)); X = X(:, shuf); b = blocks(shuf);
names = arrayfun(@(i) sprintf('v%d', i), 1:numel(b), 'UniformOutput', false);
OUT = plot_correlation_matrix(X, 'names', names, 'reorder_by_rastermap', 'nofigure');
tc.verifyEqual(sort(OUT.perm_order(:))', 1:numel(b));
tc.verifyEqual(sum(diff(b(OUT.perm_order)) ~= 0), 2);       % blocks are contiguous
% also from a precomputed r matrix
OUT2 = plot_correlation_matrix(corr(X), 'input_is_r', true, 'names', names, 'reorder_by_rastermap', 'nofigure');
tc.verifyEqual(sum(diff(b(OUT2.perm_order)) ~= 0), 2);
end
