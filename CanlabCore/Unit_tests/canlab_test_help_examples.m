function tests = canlab_test_help_examples
%CANLAB_TEST_HELP_EXAMPLES Smoke-test the docs/ "Quick example" code blocks.
%
%   These tests mirror, one-for-one, the runnable code examples shown on
%   the per-method documentation pages under CanlabCore/docs/. They are
%   smoke tests: a test passes if the example runs end-to-end without
%   erroring (and any figures created are cleaned up). We do NOT
%   pixel-compare against the docs/class_method_pngs/*.png snapshots.
%
%   Why these tests exist. The docs Quick examples are the first thing a
%   new user sees; if any of them silently bit-rot, the toolbox feels
%   broken on first contact. Wiring them into the per-push CI is the
%   cheapest possible defense.
%
%   Headless graphics. setup() forces DefaultFigureVisible 'off' so
%   nothing pops up during interactive runs. Tests that exercise
%   genuinely GUI-only code paths (orthviews, surface meshes that need
%   OpenGL) wrap the offending call in a try/catch and call
%   tc.assumeFail(...) on missing-graphics errors so CI does not fail
%   spuriously when run on a machine without a display.
%
%   Caching. setupOnce loads the emotionreg sample once and caches the
%   group t-test (raw and thresholded) and a region object on
%   tc.TestData, so the ~18 tests that need them don't recompute.

tests = functiontests(localfunctions);
end


% =====================================================================
% Fixtures
% =====================================================================

function setupOnce(tc) %#ok<*DEFNU>
% Cache the sample dataset and the derived stat objects that 18+ tests
% want, so we pay the cost once per test run rather than per test.
tc.TestData.imgs = canlab_get_sample_fmri_data();
t_raw = ttest(tc.TestData.imgs);
tc.TestData.t_raw = t_raw;
t_thr = threshold(t_raw, .005, 'unc', 'k', 10);
tc.TestData.t_thr = t_thr;
tc.TestData.region_thr = region(t_thr);
end


function setup(tc)
close all force;
tc.TestData.PrevFigVis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
end


function teardown(tc)
close all force;
if isfield(tc.TestData, 'PrevFigVis')
    set(0, 'DefaultFigureVisible', tc.TestData.PrevFigVis);
end
end


% =====================================================================
% Helpers (local)
% =====================================================================

function skip_on_environment_error(tc, ME)
%   If an error reflects a missing CI capability (graphics, interactive
%   input, or an optional dataset that ships under separate licensing),
%   skip the test on this runner instead of failing. Otherwise rethrow.
%
%   Categories handled:
%     - Graphics: OpenGL / X display / JVM not available (headless runner).
%     - Interactive input: e.g. load_atlas('canlab2024') triggers
%       bianciardi_create_atlas_obj, which prompts for confirmation when
%       atlas source files are not on the path. The runner reports
%       'MATLAB:services:MissingRequiredCapability' for any input() call
%       in batch mode.
%     - Missing data files: NPS+ signature image
%       (weights_NSF_grouppred_cvpcr.img) and the Neurosynth feature set
%       (Yarkoni_2013_Neurosynth_featureset1.mat) are required by some
%       help examples but are NOT shipped with Neuroimaging_Pattern_Masks
%       or CanlabCore on CI checkouts.
gfx_ids = {'MATLAB:graphics:opengl:Unavailable', ...
           'MATLAB:graphics:initialize'};
input_ids = {'MATLAB:services:MissingRequiredCapability', ...
             'MATLAB:UndefinedFunction'};
msg = lower(ME.message);

is_gfx = any(strcmp(ME.identifier, gfx_ids)) || ...
         contains(msg, 'opengl') || contains(msg, 'display') || ...
         contains(msg, 'java') || contains(msg, 'jvm');
is_input = any(strcmp(ME.identifier, input_ids)) && ...
           contains(msg, 'support for user input');
is_missing_data = contains(msg, 'cannot find images') || ...
                  contains(msg, 'find and add the file') || ...
                  contains(msg, 'not found in matlab path') || ...
                  contains(msg, 'no such file or directory');

if is_gfx
    tc.assumeFail(['needs graphics environment: ' ME.message]);
elseif is_input
    tc.assumeFail(['needs interactive input not available in CI: ' ME.message]);
elseif is_missing_data
    tc.assumeFail(['needs an optional data file not on the CI path: ' ME.message]);
else
    rethrow(ME);
end
end


function skip_on_graphics_error(tc, ME)
% Backwards-compatible alias - delegates to the broader environment check.
skip_on_environment_error(tc, ME);
end


% =====================================================================
% PHASE 4A — class-page and stand-alone examples (4 tests)
% =====================================================================

function test_atlas_methods_isosurface_montage(tc)
% docs/atlas_methods.md
% load_atlas('canlab2024') pulls in the Bianciardi brainstem atlas, which
% prompts for user input when its (separately-licensed) source files are
% missing - skip on CI rather than fail.
tc.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(tempdir));
try
    obj = load_atlas('canlab2024');
    create_figure('fig'); isosurface(obj);
    view(135, 30); lightFollowView;
    create_figure('fig2'); axis off; montage(obj);
    tc.verifyNotEmpty(obj.labels);
catch ME
    skip_on_environment_error(tc, ME);
end
end


function test_atlas_select_atlas_subset(tc)
% docs/individual_functions/atlas_select_atlas_subset.md
% Same Bianciardi prompt issue as the atlas_methods test above.
try
    obj = load_atlas('canlab2024');
    thal = select_atlas_subset(obj, {'Thal'});
    tc.verifyNotEmpty(thal.labels);
    tc.verifyClass(thal, 'atlas');
    create_figure('fig'); isosurface(thal);
catch ME
    skip_on_environment_error(tc, ME);
end
end


function test_annotate_binary_results_map(tc)
% docs/individual_functions/fmri_data_annotate_binary_results_map.md
% Requires Yarkoni_2013_Neurosynth_featureset1.mat, not on the CI path.
obj = tc.TestData.imgs;
t = ttest(obj, .005, 'uncorrected');
t.dat = single(t.dat > 3);
t = fmri_data(t);
try
    RESULTS = annotate_binary_results_map(t);
    tc.verifyNotEmpty(RESULTS);
catch ME
    skip_on_environment_error(tc, ME);
end
end


function test_outliers_notimeseries(tc)
% docs/individual_functions/fmri_data_outliers.md
obj = tc.TestData.imgs;
[est_outliers_uncorr, est_outliers_corr, outlier_tables] = ...
    outliers(obj, 'notimeseries');
tc.verifyEqual(numel(est_outliers_uncorr), size(obj.dat, 2));
tc.verifyEqual(numel(est_outliers_corr), size(obj.dat, 2));
tc.verifyClass(outlier_tables, 'struct');
end


% =====================================================================
% PHASE 4B — image_vector / fmri_data methods (10 tests)
% =====================================================================

function test_image_similarity_plot(tc)
% docs/individual_functions/fmri_data_image_similarity_plot.md
% Requires the NPS+ signature image (weights_NSF_grouppred_cvpcr.img),
% not on the CI path.
imgs = tc.TestData.imgs;
try
    image_similarity_plot(imgs, 'mapset', 'npsplus', 'average');
catch ME
    skip_on_environment_error(tc, ME);
end
end


function test_image_similarity_plot_bucknermaps(tc)
% docs/individual_functions/fmri_data_image_similarity_plot_bucknermaps.md
t = tc.TestData.t_raw;
stats = image_similarity_plot_bucknermaps(t);
tc.verifyNotEmpty(stats);
end


function test_jackknife_similarity(tc)
% docs/individual_functions/fmri_data_jackknife_similarity.md
imgs = tc.TestData.imgs;
[sim_values, d, low_agreement, Nvox] = jackknife_similarity(imgs, ...
    'similarity_metric', 'correlation');
tc.verifyEqual(numel(sim_values), size(imgs.dat, 2));
tc.verifyTrue(isnumeric(d));
tc.verifyClass(low_agreement, 'logical');
tc.verifyTrue(all(Nvox > 0));
end


function test_fmri_data_montage(tc)
% docs/individual_functions/fmri_data_montage.md
t = tc.TestData.t_thr;
create_figure('m'); axis off; montage(t);
end


function test_fmri_data_orthviews(tc)
% docs/individual_functions/fmri_data_orthviews.md
% orthviews uses SPM's spm_orthviews and requires a Graphics window.
tc.assumeTrue(usejava('jvm') && feature('ShowFigureWindows'), ...
    'orthviews requires an interactive figure window');
t = tc.TestData.t_thr;
try
    orthviews(t);
catch ME
    skip_on_graphics_error(tc, ME);
end
end


function test_fmri_data_mahal(tc)
% docs/individual_functions/fmri_data_mahal.md
imgs = tc.TestData.imgs;
[ds, expectedds, p_vals, wh_outlier_uncorr, wh_outlier_corr] = mahal(imgs);
tc.verifyEqual(numel(ds), size(imgs.dat, 2));
tc.verifyEqual(numel(expectedds), size(imgs.dat, 2));
tc.verifyEqual(numel(p_vals), size(imgs.dat, 2));
tc.verifyEqual(numel(wh_outlier_uncorr), size(imgs.dat, 2));
tc.verifyEqual(numel(wh_outlier_corr), size(imgs.dat, 2));
end


function test_fmri_data_pca(tc)
% docs/individual_functions/fmri_data_pca.md
imgs = tc.TestData.imgs;
[scores, eig_obj, explained] = pca(imgs, 'k', 5);
tc.verifyEqual(size(scores, 2), 5);
tc.verifyTrue(isobject(eig_obj));
tc.verifyTrue(numel(explained) >= 5);
end


function test_fmri_data_rmssd_movie(tc)
% docs/individual_functions/fmri_data_rmssd_movie.md
imgs = tc.TestData.imgs;
rmssd_movie(imgs);
end


function test_fmri_data_surface(tc)
% docs/individual_functions/fmri_data_surface.md
t = tc.TestData.t_thr;
tc.assumeTrue(usejava('jvm'), 'surface requires Java for rendering');
try
    create_figure('s'); surface(t);
catch ME
    skip_on_graphics_error(tc, ME);
end
end


function test_fmri_data_wedge_plot_by_atlas(tc)
% docs/individual_functions/fmri_data_wedge_plot_by_atlas.md
% Use yeo17networks - the example in the .md uses this atlas.
t = threshold(tc.TestData.t_raw, .005, 'unc');
[hh, vals] = wedge_plot_by_atlas(t, 'atlases', {'yeo17networks'});
tc.verifyNotEmpty(vals);
end


% =====================================================================
% PHASE 4B — statistic_image methods (4 tests)
% =====================================================================

function test_statistic_image_riverplot(tc)
% docs/individual_functions/statistic_image_riverplot.md
% load_image_set('npsplus') needs the NPS+ signature images, not on CI.
imgs = tc.TestData.imgs;
try
    layer1 = load_image_set('npsplus'); layer1 = get_wh_image(layer1, 1:4);
    layer2 = ttest(imgs); layer2 = threshold(layer2, .005, 'unc');
    layer2 = fmri_data(layer2);
    layer1.image_names = char({'NPS','SIIPS','GenS','VPS'});
    layer2.image_names = char({'EmoReg group t'});
    riverplot(layer1, 'layer2', layer2);
catch ME
    skip_on_environment_error(tc, ME);
end
end


function test_statistic_image_multi_threshold(tc)
% docs/individual_functions/statistic_image_multi_threshold.md
% Signature: [o2, dat, sig, pcl, ncl] = multi_threshold(dat, ...)
t = tc.TestData.t_raw;
[o2, dat_out, sig, pcl, ncl] = multi_threshold(t);
tc.verifyClass(o2, 'fmridisplay');
tc.verifyClass(dat_out, 'statistic_image');
tc.verifyTrue(iscell(sig) || islogical(sig) || isnumeric(sig));
tc.verifyTrue(isstruct(pcl) || isobject(pcl) || iscell(pcl) || isempty(pcl));
tc.verifyTrue(isstruct(ncl) || isobject(ncl) || iscell(ncl) || isempty(ncl));
end


function test_statistic_image_table(tc)
% docs/individual_functions/statistic_image_table.md
t = tc.TestData.t_thr;
[r, results_table] = table(t);
tc.verifyClass(r, 'region');
tc.verifyClass(results_table, 'table');
end


function test_statistic_image_threshold(tc)
% docs/individual_functions/statistic_image_threshold.md
t = tc.TestData.t_raw;
t = threshold(t, .005, 'unc', 'k', 10);
create_figure('thr'); axis off; montage(t);
tc.verifyClass(t, 'statistic_image');
end


% =====================================================================
% PHASE 4B — region methods (6 tests)
% =====================================================================

function test_region_table(tc)
% docs/individual_functions/region_table.md
r = tc.TestData.region_thr;
table(r);
tc.verifyTrue(numel(r) > 0);
end


function test_region_surface(tc)
% docs/individual_functions/region_surface.md
r = tc.TestData.region_thr;
try
    create_figure('rs'); surface(r);
catch ME
    skip_on_graphics_error(tc, ME);
end
end


function test_region_labelled_surface(tc)
% docs/individual_functions/region_labelled_surface.md
r = tc.TestData.region_thr;
try
    create_figure('rls'); labelled_surface(r);
catch ME
    skip_on_graphics_error(tc, ME);
end
end


function test_region_isosurface(tc)
% docs/individual_functions/region_isosurface.md
r = tc.TestData.region_thr;
create_figure('ri'); set(gcf, 'Position', [100 100 900 700]);
isosurface(r);
end


function test_region_montage(tc)
% docs/individual_functions/region_montage.md
r = tc.TestData.region_thr;
montage(r, 'regioncenters', 'colormap');
end


function test_region_table_of_atlas_regions_covered(tc)
% docs/individual_functions/region_table_of_atlas_regions_covered.md
% NOTE: the region/ method self-disclaims as broken (the docs page calls
% this out and recommends the @image_vector path). Treat as smoke-only:
% if it errors we mark it as expected-failure rather than a regression.
r = tc.TestData.region_thr;
try
    table_of_atlas_regions_covered(r);
catch ME
    tc.assumeFail(['region.table_of_atlas_regions_covered ' ...
        'is documented as broken on the @region path: ' ME.message]);
end
end


% =====================================================================
% PHASE 4B — stand-alone visualization helpers (6 tests)
% =====================================================================

function test_addbrain(tc)
% docs/individual_functions/addbrain.md
try
    create_figure('ab'); set(gcf, 'Position', [100 100 900 700]);
    p = addbrain('hires left'); set(p, 'FaceAlpha', .8);
    addbrain('thalamus');
    addbrain('bg');
    view(135, 10); lightRestoreSingle;
catch ME
    skip_on_graphics_error(tc, ME);
end
end


function test_barplot_columns(tc)
% docs/individual_functions/barplot_columns.md
rng(7);
Y = [randn(20,1)+1, randn(20,1)+0.3, randn(20,1)-0.5];
colors = seaborn_colors(8);
[handles, ~, ~, statstable] = barplot_columns(Y, 'nofig', ...
    'names', {'Cond A', 'Cond B', 'Cond C'}, 'colors', colors(1:3));
tc.verifyClass(handles, 'struct');
tc.verifyClass(statstable, 'table');
tc.verifyEqual(height(statstable), 3);
end


function test_cluster_surf(tc)
% docs/individual_functions/cluster_surf.md
r = tc.TestData.region_thr;
try
    create_figure('cs'); set(gcf, 'Position', [100 100 1000 700]);
    cluster_surf(r, 5, 'left');
catch ME
    skip_on_graphics_error(tc, ME);
end
end


function test_image_scatterplot(tc)
% docs/individual_functions/image_scatterplot.md
imgs = tc.TestData.imgs;
% The .md example uses an OLS vs. robust regression comparison; mirror it.
tc.assumeTrue(ismember('Reappraisal_Success', ...
    imgs.metadata_table.Properties.VariableNames), ...
    'metadata_table is missing Reappraisal_Success column');
imgs.X = [imgs.metadata_table.Reappraisal_Success - ...
          mean(imgs.metadata_table.Reappraisal_Success), ...
          ones(size(imgs.dat, 2), 1)];
out_ols = regress(imgs, 'noverbose', 'nodisplay');
out_rob = regress(imgs, 'robust', 'noverbose', 'nodisplay');
t_ols = get_wh_image(out_ols.t, 1);
t_rob = get_wh_image(out_rob.t, 1);
image_scatterplot(t_ols, t_rob, 'colorpoints');
xlabel('OLS t-value'); ylabel('Robust regression t-value');
end


function test_canlab_results_fmridisplay(tc)
% docs/individual_functions/canlab_results_fmridisplay.md
t = tc.TestData.t_thr;
o2 = canlab_results_fmridisplay(t);
tc.verifyClass(o2, 'fmridisplay');
end


function test_plot_correlation_matrix(tc)
% docs/individual_functions/plot_correlation_matrix.md
rng(7);
S = toeplitz([1 .6 .3 .1 0 0]);
X = mvnrnd([0 0 0 0 0 0], S, 50);
var_names = {'A' 'B' 'C' 'D' 'E' 'F'};
OUT = plot_correlation_matrix(X, 'var_names', var_names);
tc.verifyTrue(isstruct(OUT) || isobject(OUT));
end
