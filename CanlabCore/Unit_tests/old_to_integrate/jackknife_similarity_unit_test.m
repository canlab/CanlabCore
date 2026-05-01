%% test_jackknife_similarity_documented.m
% TEST_JACKKNIFE_SIMILARITY_DOCUMENTED
%
% This script provides a more systematic and interpretable set of tests for
% jackknife_similarity.m using the sample image set returned by:
%
%     obj = load_image_set('emotionreg');
%
% The loaded object is an fmri_data object containing 30 images. The goals
% of the tests are:
%
%   1. Unit-test the main options and metrics to make sure outputs have the
%      expected size/type and run without errors on a standard dataset.
%   2. Compare how different image rescaling / normalization procedures
%      affect the mean agreement value and agreement effect size (d) across
%      all metrics implemented in jackknife_similarity.
%   3. Compare how additive image-wise location shifts affect agreement for
%      all metrics.
%   4. Compare how multiplicative image-wise scale shifts affect agreement
%      for all metrics.
%
% Rationale
% ---------
% Different agreement metrics answer different questions:
%
% - Correlation primarily captures spatial pattern similarity after
%   centering and scaling within each image pair.
% - Cosine similarity is sensitive to mean shifts because it does not
%   center the vectors.
% - Concordance correlation and absolute agreement are sensitive to
%   absolute value differences, including location and scale changes.
% - standardized_abs_deviation, mean_shift_z, and scale_shift_z were added
%   specifically to capture image-wise deviations in absolute level and
%   scale relative to the other images.
%
% This script therefore treats location and scale perturbations as
% meaningful structured error and asks which metrics are more or less
% sensitive to them.
%
% Notes
% -----
% - The script uses only jackknife_similarity outputs:
%       sim_values, Nvox, d, low_agreement
% - All comparisons are summarized in tables and interpreted at the end.
% - Random perturbations are reproducible using a fixed RNG seed.
%
% -------------------------------------------------------------------------

clear; clc;
rng(1, 'twister');

%% Setup
% Load the standard CANlab sample image set. This is the baseline dataset
% used throughout the script.

obj = load_image_set('emotionreg');
n_images = size(obj.dat, 2);

fprintf('Loaded emotionreg image set with %d images.\n', n_images);

% Metrics to compare. These include the newer deviation-based metrics.
metrics = { ...
    'correlation' ...
    'cosine_similarity' ...
    'dot_product' ...
    'normalized_absolute_agreement' ...
    'concordance_correlation' ...
    'standardized_abs_deviation' ...
    'mean_shift_z' ...
    'scale_shift_z'};

n_metrics = numel(metrics);

%% 1) Unit tests on metrics and key options
% Rationale:
% This section verifies that each metric runs on the sample dataset and
% returns outputs with expected dimensions. It also checks a few key option
% combinations that can change behavior:
%
%   - default options
%   - complete_cases = true
%   - treat_zero_as_data = true
%   - verbose = false
%   - plot keyword
%
% The goal is not to prove numerical correctness in a formal software-
% engineering sense, but to catch obvious failures and ensure the interface
% works as expected on a representative fmri_data object.

fprintf('\n============================================================\n');
fprintf('SECTION 1: Unit tests of metrics and key options\n');
fprintf('============================================================\n');

unit_results = table('Size', [n_metrics 7], ...
    'VariableTypes', {'string','double','double','double','double','double','double'}, ...
    'VariableNames', {'metric','mean_similarity','d','n_low_agreement', ...
                      'mean_Nvox','min_similarity','max_similarity'});

for m = 1:n_metrics
    metric = metrics{m};

    [s, d, low_agreement, Nvox] = jackknife_similarity(obj, ...
        'similarity_metric', metric, ...
        'verbose', false);

    assert(isvector(s) && numel(s) == n_images, ...
        'Similarity output has wrong size for metric: %s', metric);
    assert(isvector(Nvox) && numel(Nvox) == n_images, ...
        'Nvox output has wrong size for metric: %s', metric);
    assert(isscalar(d), ...
        'd output is not scalar for metric: %s', metric);
    assert(islogical(low_agreement) && numel(low_agreement) == n_images, ...
        'low_agreement output has wrong size/type for metric: %s', metric);

    % Additional option smoke tests
    [~, ~, ~, ~] = jackknife_similarity(obj, ...
        'similarity_metric', metric, ...
        'complete_cases', true, ...
        'verbose', false);

    [~, ~, ~, ~] = jackknife_similarity(obj, ...
        'similarity_metric', metric, ...
        'treat_zero_as_data', true, ...
        'verbose', false);

    [~, ~, ~, ~] = jackknife_similarity(obj, ...
        'similarity_metric', metric, ...
        'verbose', false);

    unit_results.metric(m) = string(metric);
    unit_results.mean_similarity(m) = mean(s, 'omitnan');
    unit_results.d(m) = d;
    unit_results.n_low_agreement(m) = sum(low_agreement);
    unit_results.mean_Nvox(m) = mean(Nvox, 'omitnan');
    unit_results.min_similarity(m) = min(s, [], 'omitnan');
    unit_results.max_similarity(m) = max(s, [], 'omitnan');
end

% One simple smoke test of the 'plot' keyword on the default metric.
% This should create a figure without error.
close all;
jackknife_similarity(obj, 'plot', 'verbose', false);
drawnow;

disp(unit_results);

%% 2) Effect of rescaling / normalization on agreement metrics
% Rationale:
% Agreement metrics may respond very differently depending on how images are
% normalized beforehand. This matters because common preprocessing choices
% can intentionally remove or emphasize:
%
%   - global mean shifts across images
%   - global scale differences across images
%   - overall tissue-intensity normalization
%
% Here we compare several transformations:
%
%   raw
%       Original images.
%
%   gm_wm_csf_norm
%       CANlab normalization using normalize_gm_by_wm_csf, if available.
%
%   zscoreimages
%       CANlab image-wise z-scoring, which removes mean and scales each
%       image to unit variance.
%
%   center_each_image
%       Removes each image's mean only.
%
%   center_and_madscale_each_image
%       Removes each image's mean and scales by each image's MAD. This is a
%       robust within-image normalization.
%
% These manipulations test which metrics primarily reflect pattern shape,
% absolute magnitude, location shifts, or scale shifts.

fprintf('\n============================================================\n');
fprintf('SECTION 2: Effect of image normalization / rescaling\n');
fprintf('============================================================\n');

datasets = struct();
datasets.raw = obj;

datasets.gm_wm_csf_norm = normalize_gm_by_wm_csf(obj);
datasets.zscoreimages = rescale(obj, 'zscoreimages');

datasets.center_each_image = obj;
datasets.center_each_image.dat = center_columns(obj.dat);

datasets.center_and_madscale_each_image = obj;
datasets.center_and_madscale_each_image.dat = robust_center_and_scale_columns(obj.dat);

dataset_names = fieldnames(datasets);
n_datasets = numel(dataset_names);

norm_results = table();
row = 0;

for di = 1:n_datasets
    ds_name = dataset_names{di};
    ds_obj = datasets.(ds_name);

    for m = 1:n_metrics
        metric = metrics{m};

        [s, d, low_agreement, Nvox] = jackknife_similarity(ds_obj, ...
            'similarity_metric', metric, ...
            'verbose', false);

        row = row + 1;
        norm_results.dataset(row,1) = string(ds_name);
        norm_results.metric(row,1) = string(metric);
        norm_results.mean_similarity(row,1) = mean(s, 'omitnan');
        norm_results.d(row,1) = d;
        norm_results.mean_Nvox(row,1) = mean(Nvox, 'omitnan');
        norm_results.n_low_agreement(row,1) = sum(low_agreement);
    end
end

disp(norm_results);

% Summaries: highest and lowest mean agreement and d for each metric
norm_summary = summarize_best_worst(norm_results, 'dataset');
disp(norm_summary);

%% 3) Effect of global additive shift (all images)
% Rationale:
% Adds a constant to ALL images equally. This should NOT affect correlation,
% but SHOULD affect metrics sensitive to absolute magnitude.

fprintf('\n============================================================\n');
fprintf('SECTION 3: Global additive shift (all images)\n');
fprintf('============================================================\n');

c = .2 * std(obj.dat(:),'omitnan');   % small global shift

obj_shift_global = obj;
obj_shift_global.dat = obj.dat + c;

results_global_shift = compare_two_datasets(obj, obj_shift_global, metrics, 'raw', 'global_shift');
disp(results_global_shift);

% Collect full similarity values
sDiff = [];

for m = 1:numel(metrics)
    metric = metrics{m};

    [s1, ~] = jackknife_similarity(obj, 'similarity_metric', metric, 'verbose', false);
    [s2, ~] = jackknife_similarity(obj_shift_global, 'similarity_metric', metric, 'verbose', false);

    sDiff(:,m) = s2 - s1;
end

create_figure('global_shift_violin',1,1);
barplot_columns(sDiff, 'colors', seaborn_colors(numel(metrics)), 'nobars')
ylabel('<--Lower agreement   Shift effect   Higher agreement-->')
set(gca,'XTickLabel',metrics,'XTickLabelRotation',45)
yline(0,'--','No change')

%% 4) Effect of global multiplicative scale (all images)
% Rationale:
% Multiply all images by same factor. Should not affect correlation, but
% should affect magnitude-sensitive metrics.

fprintf('\n============================================================\n');
fprintf('SECTION 4: Global scale (all images)\n');
fprintf('============================================================\n');

a = 1.5;

obj_scale_global = obj;
obj_scale_global.dat = obj.dat * a;

results_global_scale = compare_two_datasets(obj, obj_scale_global, metrics, 'raw', 'global_scale');
disp(results_global_scale);

% Collect full similarity values
sDiff = [];

for m = 1:numel(metrics)
    metric = metrics{m};

    [s1, ~] = jackknife_similarity(obj, 'similarity_metric', metric, 'verbose', false);
    [s2, ~] = jackknife_similarity(obj_scale_global, 'similarity_metric', metric, 'verbose', false);

    sDiff(:,m) = s2 - s1;
end

create_figure('global_scale_violin',1,1);
barplot_columns(sDiff, 'colors', seaborn_colors(numel(metrics)), 'nobars')
ylabel('<--Lower agreement   Scale effect   Higher agreement-->')
set(gca,'XTickLabel',metrics,'XTickLabelRotation',45)
yline(0,'--','No change')

%% 5) Effect of additive image-wise location shifts
% Rationale:
% This section introduces random image-specific additive offsets:
%
%     X_i -> X_i + c_i
%
% where c_i differs across images but is constant within each image.
%
% This models structured error in image location / baseline intensity.
% Since the offsets differ across images, these are not expected to be
% removed by all metrics. We ask:
%
%   - Which metrics are more sensitive to additive offsets?
%   - Do mean agreement and d increase or decrease?
%
% The shift magnitude is tied to the empirical image-wise SD of image means
% so that the perturbation is noticeable but not pathological.

fprintf('\n============================================================\n');
fprintf('SECTION 5: Effect of additive image-wise location shifts\n');
fprintf('============================================================\n');

obj_shift = obj;
img_means = mean(obj.dat, 1, 'omitnan');
shift_sigma = std(img_means, 0, 2, 'omitnan');

if ~isfinite(shift_sigma) || shift_sigma == 0
    shift_sigma = median(abs(obj.dat(~isnan(obj.dat))));
end

location_offsets = randn(1, n_images) .* shift_sigma;
obj_shift.dat = obj_shift.dat + location_offsets;

location_results = compare_two_datasets(obj, obj_shift, metrics, 'raw', 'location_shift');
disp(location_results);

location_summary = summarize_shift_effects(location_results);
disp(location_summary);

% Full similarity differences
sDiff = [];

for m = 1:numel(metrics)
    metric = metrics{m};

    [s1, ~] = jackknife_similarity(obj, 'similarity_metric', metric, 'verbose', false);
    [s2, ~] = jackknife_similarity(obj_shift, 'similarity_metric', metric, 'verbose', false);

    sDiff(:,m) = s2 - s1;
end

create_figure('random_shift_violin',1,1);
barplot_columns(sDiff, 'colors', seaborn_colors(numel(metrics)), 'nobars')
ylabel('<--Lower agreement   Random shift   Higher agreement-->')
set(gca,'XTickLabel',metrics,'XTickLabelRotation',45)
yline(0,'--','No change')

%% 6) Effect of multiplicative image-wise scale shifts
% Rationale:
% This section introduces random positive multiplicative coefficients:
%
%     X_i -> a_i X_i,   with a_i > 0
%
% This models structured scale variability across images. It preserves the
% sign pattern of each image but changes the overall magnitude image by
% image. We again compare mean agreement and d across metrics.
%
% The coefficients are lognormally distributed so that all multipliers are
% positive and the perturbation is approximately symmetric on the log scale.

fprintf('\n============================================================\n');
fprintf('SECTION 6: Effect of multiplicative image-wise scale shifts\n');
fprintf('============================================================\n');

obj_scale = obj;
scale_sigma = 0.5;                        % moderate multiplicative spread
obj_scale = add_scale_noise_preserve_mean(obj, scale_sigma);

% scale_factors = exp(scale_sigma .* randn(1, n_images));
% obj_scale.dat = obj_scale.dat .* scale_factors;

scale_results = compare_two_datasets(obj, obj_scale, metrics, 'raw', 'scale_shift');
disp(scale_results);

scale_summary = summarize_shift_effects(scale_results);
disp(scale_summary);

sDiff = [];

for m = 1:numel(metrics)
    metric = metrics{m};

    [s1, ~] = jackknife_similarity(obj, 'similarity_metric', metric, 'verbose', false);
    [s2, ~] = jackknife_similarity(obj_scale, 'similarity_metric', metric, 'verbose', false);

    sDiff(:,m) = s2 - s1;
end

create_figure('random_scale_violin',1,1);
barplot_columns(sDiff, 'colors', seaborn_colors(numel(metrics)), 'nobars')
ylabel('<--Lower agreement   Random scale   Higher agreement-->')
set(gca,'XTickLabel',metrics,'XTickLabelRotation',45)
yline(0,'--','No change')

%% 7) Effect of voxelwise Gaussian noise
% Rationale:
% Adds independent voxel noise. Should reduce ALL agreement metrics.

fprintf('\n============================================================\n');
fprintf('SECTION 7: Pattern noise (voxelwise Gaussian)\n');
fprintf('============================================================\n');

noise_sigma = 1.5 * std(obj.dat(:),'omitnan');
obj_noise = add_pattern_noise_preserve_scale(obj, noise_sigma);

% obj_noise = obj;
% obj_noise.dat = obj.dat + noise_sigma * randn(size(obj.dat));

results_noise = compare_two_datasets(obj, obj_noise, metrics, 'raw', 'noise');
disp(results_noise);

% Collect similarity differences
sDiff = [];

for m = 1:numel(metrics)
    metric = metrics{m};

    [s1, ~] = jackknife_similarity(obj, 'similarity_metric', metric, 'verbose', false);
    [s2, ~] = jackknife_similarity(obj_noise, 'similarity_metric', metric, 'verbose', false);

    sDiff(:,m) = s2 - s1;
end

create_figure('noise_violin',1,1);
barplot_columns(sDiff, 'colors', seaborn_colors(numel(metrics)), 'nobars')
ylabel('<--Lower agreement   Noise effect   Higher agreement-->')
set(gca,'XTickLabel',metrics,'XTickLabelRotation',45)
yline(0,'--','No change')

%% FINAL SUMMARY: Heatmap of sensitivity

fprintf('\n============================================================\n');
fprintf('FINAL SUMMARY: Metric sensitivity heatmap\n');
fprintf('============================================================\n');

redblue = colormap_tor([0 0 1], [1 0 0], [1 1 1]);

% Collect delta mean + delta d across all manipulations

all_results = { ...
    results_global_shift ...
    results_global_scale ...
    location_results ...
    scale_results ...
    results_noise };

labels = {'global shift','global scale','random shift','random scale','pattern noise'};

n_tests = numel(all_results);
n_metrics = numel(metrics);

delta_mean = zeros(n_metrics, n_tests);
delta_d = zeros(n_metrics, n_tests);

for t = 1:n_tests
    tbl = all_results{t};

    for m = 1:n_metrics
        delta_mean(m,t) = tbl.delta_mean_similarity(m);
        delta_d(m,t) = tbl.delta_d(m);
    end
end

% Normalize to remove differences in scale across metrics and noise
% manipulations, keeping the zero-point intact (we expect decreases not
% increases)

Mnorm = delta_mean;

for i = 1:5
    Mnorm = Mnorm ./ max(abs(Mnorm), [], 2);  % normalize rows
    Mnorm = Mnorm ./ max(abs(Mnorm), [], 1);  % normalize columns
end

Dnorm = delta_d;

for i = 1:5
    Dnorm = Dnorm ./ max(abs(Dnorm), [], 2);  % normalize rows
    Dnorm = Dnorm ./ max(abs(Dnorm), [], 1);  % normalize columns
end

% Eliminate very small effects for clarity
Dnorm(abs(delta_d) < 0.2) = 0;

create_figure('sensitivity_heatmap');

imagesc(Mnorm)
title('Mean agreement sensitivity')
set(gca,'XTick',1:n_tests,'XTickLabel',labels,...
    'YTick',1:n_metrics,'YTickLabel', format_strings_for_legend(metrics))
colorbar
colormap(redblue)   % CANlab colormap
set(gca, 'CLim', [-max(abs(get(gca, 'CLim'))) max(abs(get(gca, 'CLim')))])
set(gca, 'YDir', 'reverse'); axis tight
xlabel('Noise Type'), ylabel('Agreement metric')
set(gcf, 'Position', [30   294   684   701])

create_figure('sensitivity_heatmap2');

% subplot(1,2,1)
imagesc(delta_d)
title('Effect size (d) change, blue=exp dir')
set(gca,'XTick',1:n_tests,'XTickLabel',labels,...
    'YTick',1:n_metrics,'YTickLabel', format_strings_for_legend(metrics))
colorbar
colormap(redblue)
set(gca, 'CLim', [-max(abs(get(gca, 'CLim'))) max(abs(get(gca, 'CLim')))])
set(gca, 'YDir', 'reverse'); axis tight
xlabel('Noise Type'), ylabel('Agreement metric')
set(gcf, 'Position', [30+500   294   684   701])

create_figure('sensitivity_heatmap3');

% subplot(1,2,2)
imagesc(Dnorm)
title('D double-normed change, blue=exp dir')
set(gca,'XTick',1:n_tests,'XTickLabel',labels,...
    'YTick',1:n_metrics,'YTickLabel', format_strings_for_legend(metrics))
colorbar
colormap(redblue)
set(gca, 'CLim', [-max(abs(get(gca, 'CLim'))) max(abs(get(gca, 'CLim')))])
set(gca, 'YDir', 'reverse'); axis tight
xlabel('Noise Type'), ylabel('Agreement metric')
set(gcf, 'Position', [30+1000   294   684   701])
%% Optional summary plots
% Compact plots to visualize the mean-agreement and d changes produced by
% the location and scale manipulations.
% 
% create_figure('jackknife_similarity_tests', 2, 2);
% 
% subplot(2,2,1);
% plot_metric_changes(location_results, 'delta_mean_similarity', 'Location shift: change in mean agreement');
% 
% subplot(2,2,2);
% plot_metric_changes(location_results, 'delta_d', 'Location shift: change in d');
% 
% subplot(2,2,3);
% plot_metric_changes(scale_results, 'delta_mean_similarity', 'Scale shift: change in mean agreement');
% 
% subplot(2,2,4);
% plot_metric_changes(scale_results, 'delta_d', 'Scale shift: change in d');

%% Interpretation of findings
% The interpretation below is generated from the observed results, not from
% hard-coded expectations. This makes the script a self-contained analysis
% report rather than just a collection of tests.

fprintf('\n============================================================\n');
fprintf('INTERPRETATION OF FINDINGS\n');
fprintf('============================================================\n');

% Section 1
fprintf('\n1) Unit tests\n');
fprintf(['All listed metrics were run on the 30-image emotionreg dataset and\n' ...
         'checked for basic output validity (vector sizes, scalar d, and\n' ...
         'logical low_agreement flags). The ''plot'' keyword was also smoke-\n' ...
         'tested. Review unit_results above for the baseline scale/range of\n' ...
         'each metric on the sample dataset.\n']);

% Section 2
fprintf('\n2) Rescaling / normalization effects\n');
for m = 1:n_metrics
    metric = metrics{m};
    wh = strcmp(norm_results.metric, metric);
    this = norm_results(wh, :);

    [mx, imx] = max(this.mean_similarity);
    [mn, imn] = min(this.mean_similarity);
    [dx, idx] = max(this.d);
    [dn, idn] = min(this.d);

    fprintf(['Metric %-28s highest mean agreement: %-28s (%8.4f), ' ...
             'lowest mean agreement: %-28s (%8.4f)\n'], ...
        metric, this.dataset(imx), mx, this.dataset(imn), mn);

    fprintf(['%31s highest d: %-28s (%8.4f), ' ...
             'lowest d: %-28s (%8.4f)\n'], ...
        '', this.dataset(idx), dx, this.dataset(idn), dn);
end

fprintf(['Interpretation: metrics whose mean agreement and d change strongly\n' ...
         'across normalization schemes are more sensitive to image location\n' ...
         'and scale. Metrics that change little are more dominated by shared\n' ...
         'spatial pattern after within-image normalization.\n']);

% Section 3
fprintf('\n3) Additive location shifts\n');
report_shift_interpretation(location_summary, 'location shifts');

fprintf(['Interpretation: if a metric decreases strongly after random image-\n' ...
         'wise offsets are added, it is sensitive to location shifts. If it\n' ...
         'changes little, it is relatively insensitive to this type of error.\n']);

% Section 4
fprintf('\n4) Multiplicative scale shifts\n');
report_shift_interpretation(scale_summary, 'scale shifts');

fprintf(['Interpretation: if a metric decreases strongly after random image-\n' ...
         'wise positive rescaling, it is sensitive to multiplicative scale\n' ...
         'differences. Metrics designed to detect scale mismatch should show\n' ...
         'larger changes in mean agreement and/or d.\n']);

fprintf('\nDone.\n');

%% Local helper functions

function Xc = center_columns(X)
% Remove column means, omitting NaNs
col_mean = mean(X, 1, 'omitnan');
Xc = X - col_mean;
end

function Xs = robust_center_and_scale_columns(X)
% Remove column means and divide each column by its MAD across voxels.
col_mean = mean(X, 1, 'omitnan');
Xc = X - col_mean;
col_mad = median(abs(Xc), 1, 'omitnan');
col_mad(col_mad == 0 | isnan(col_mad)) = 1;
Xs = Xc ./ col_mad;
end

function out = compare_two_datasets(obj1, obj2, metrics, label1, label2)
% Compare two datasets metric-by-metric using mean agreement and d.
n_metrics = numel(metrics);
out = table('Size', [n_metrics 7], ...
    'VariableTypes', {'string','double','double','double','double','double','double'}, ...
    'VariableNames', {'metric','mean_similarity_1','mean_similarity_2', ...
                      'd_1','d_2','delta_mean_similarity','delta_d'});

for m = 1:n_metrics
    metric = metrics{m};

    % [s, d, low_agreement, Nvox]
    [s1, d1] = jackknife_similarity(obj1, 'similarity_metric', metric, 'verbose', false);
    [s2, d2] = jackknife_similarity(obj2, 'similarity_metric', metric, 'verbose', false);

    out.metric(m) = string(metric);
    out.mean_similarity_1(m) = mean(s1, 'omitnan');
    out.mean_similarity_2(m) = mean(s2, 'omitnan');
    out.d_1(m) = d1;
    out.d_2(m) = d2;
    out.delta_mean_similarity(m) = out.mean_similarity_2(m) - out.mean_similarity_1(m);
    out.delta_d(m) = out.d_2(m) - out.d_1(m);
end

out.Properties.Description = sprintf('Comparison: %s vs %s', label1, label2);
end

function out = summarize_best_worst(tbl, groupvar)
% For each metric, find highest/lowest mean agreement and d.
metrics = unique(tbl.metric);
out = table();

for i = 1:numel(metrics)
    wh = strcmp(tbl.metric, metrics(i));
    t = tbl(wh, :);

    [mx, imx] = max(t.mean_similarity);
    [mn, imn] = min(t.mean_similarity);
    [dx, idx] = max(t.d);
    [dn, idn] = min(t.d);

    out.metric(i,1) = metrics(i);
    out.highest_mean_source(i,1) = t.(groupvar)(imx);
    out.highest_mean_similarity(i,1) = mx;
    out.lowest_mean_source(i,1) = t.(groupvar)(imn);
    out.lowest_mean_similarity(i,1) = mn;

    out.highest_d_source(i,1) = t.(groupvar)(idx);
    out.highest_d(i,1) = dx;
    out.lowest_d_source(i,1) = t.(groupvar)(idn);
    out.lowest_d(i,1) = dn;
end
end

function out = summarize_shift_effects(tbl)
% Rank metrics by how much mean agreement and d changed.
[~, ix_mean_down] = sort(tbl.delta_mean_similarity, 'ascend');
[~, ix_mean_up] = sort(tbl.delta_mean_similarity, 'descend');
[~, ix_d_down] = sort(tbl.delta_d, 'ascend');
[~, ix_d_up] = sort(tbl.delta_d, 'descend');

out = struct();
out.most_decreased_mean_metric = tbl.metric(ix_mean_down(1));
out.most_decreased_mean_value = tbl.delta_mean_similarity(ix_mean_down(1));
out.most_increased_mean_metric = tbl.metric(ix_mean_up(1));
out.most_increased_mean_value = tbl.delta_mean_similarity(ix_mean_up(1));

out.most_decreased_d_metric = tbl.metric(ix_d_down(1));
out.most_decreased_d_value = tbl.delta_d(ix_d_down(1));
out.most_increased_d_metric = tbl.metric(ix_d_up(1));
out.most_increased_d_value = tbl.delta_d(ix_d_up(1));
end

function plot_metric_changes(tbl, varname, ttl)
vals = tbl.(varname);
bar(vals);
set(gca, 'XTick', 1:numel(tbl.metric), 'XTickLabel', cellstr(tbl.metric), ...
    'XTickLabelRotation', 45);
ylabel(varname, 'Interpreter', 'none');
title(ttl, 'Interpreter', 'none');
yline(0, '--');
box on
end

function report_shift_interpretation(summary_struct, label)

fprintf('Largest decrease in mean agreement under %s: %s (%8.4f)\n', ...
    label, summary_struct.most_decreased_mean_metric, summary_struct.most_decreased_mean_value);

fprintf('Largest increase in mean agreement under %s: %s (%8.4f)\n', ...
    label, summary_struct.most_increased_mean_metric, summary_struct.most_increased_mean_value);

fprintf('Largest decrease in d under %s: %s (%8.4f)\n', ...
    label, summary_struct.most_decreased_d_metric, summary_struct.most_decreased_d_value);

fprintf('Largest increase in d under %s: %s (%8.4f)\n', ...
    label, summary_struct.most_increased_d_metric, summary_struct.most_increased_d_value);

end

% Just multiplying by a random value changes the mean of the image too.
% We want to preserve the mean here and just alter the scale.

function obj_scale = add_scale_noise_preserve_mean(obj, scale_sigma)

dat = obj.dat;
[nvox, nimg] = size(dat);

obj_scale = obj;

for i = 1:nimg

    x = dat(:, i);

    mu = mean(x, 'omitnan');

    % center
    xc = x - mu;

    % draw multiplicative factor (lognormal is nice)
    a = exp(scale_sigma * randn);

    % scale about mean
    x_scaled = xc * a + mu;

    obj_scale.dat(:, i) = x_scaled;

end

end % function

% Just adding Gaussian noise changes the scale of the image too, increasing
% the scale. We want to preserve the scale here and just alter the pattern.

function obj_noise = add_pattern_noise_preserve_scale(obj, noise_level)

dat = obj.dat;
[nvox, nimg] = size(dat);

obj_noise = obj;

for i = 1:nimg

    x = dat(:, i);

    % original stats
    mu = mean(x, 'omitnan');
    mad0 = median(abs(x - mu), 'omitnan');

    if mad0 == 0 || ~isfinite(mad0)
        continue
    end

    % add Gaussian noise (scaled to original MAD)
    noise = noise_level * mad0 * randn(size(x));
    x_noisy = x + noise;

    % re-center
    mu_new = mean(x_noisy, 'omitnan');
    x_noisy = x_noisy - mu_new + mu;

    % rescale to original MAD
    mad_new = median(abs(x_noisy - mu), 'omitnan');

    if mad_new > 0
        x_noisy = (x_noisy - mu) * (mad0 / mad_new) + mu;
    end

    obj_noise.dat(:, i) = x_noisy;

end

end % function

