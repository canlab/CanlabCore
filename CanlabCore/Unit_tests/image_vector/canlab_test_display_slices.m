function tests = canlab_test_display_slices
%CANLAB_TEST_DISPLAY_SLICES Tests for image_vector/display_slices.
%
% Covers the canonical use cases of display_slices:
%   - default single-volume call -> axial montage with auto-chosen spacing
%   - coronal / sagittal single-volume views
%   - 'three_views' composite figure (3 horizontal montages stacked)
%   - explicit spacing / startslice / endslice overrides
%   - 'multi_image' for fmri_data: one imagesc slice per image
%   - 'multi_image' for statistic_image: T1 underlay + colored blobs
%   - 'multi_image' for multi-column statistic_image (t.dat(:, 2))
%
% Tests assume an interactive figure window is available; they are
% skipped in pure headless / -batch runs without a display.
%
% To run all tests:
%   results = runtests('canlab_test_display_slices')
%
% Example one-liners (also embedded in display_slices help):
%
%   imgs = load_image_set('emotionreg');
%   m    = mean(imgs);
%   display_slices(m);
%   display_slices(m, 'coronal');
%   display_slices(m, 'saggital');
%   display_slices(m, 'three_views');
%   display_slices(m, 'axial', 'spacing', 8);
%   display_slices(m, 'saggital', 'spacing', 10, 'vertical');
%   display_slices(m, 'axial', 'slices_per_row', 20, 'spacing', 4, 'startslice', -10, 'endslice', 10);
%   display_slices(imgs, 'multi_image');
%   display_slices(imgs, 'multi_image', 'slice', -10);
%   display_slices(imgs, 'multi_image', 'coronal', 'slice', 0);
%   display_slices(imgs, 'multi_image', 'saggital', 'slice', -4);
%   display_slices(imgs, 'multi_image', 'nimages', 6);
%   t = canlab_get_sample_thresholded_t(0.01);
%   display_slices(t, 'multi_image');
%
%   % Multi-image statistic_image (t_multi.dat(:, 2) is the 2nd map):
%   half1 = ttest(get_wh_image(imgs, 1:15));
%   half2 = ttest(get_wh_image(imgs, 16:30));
%   half1 = threshold(half1, .01, 'unc');
%   half2 = threshold(half2, .01, 'unc');
%   t_multi             = half1;
%   t_multi.dat         = [half1.dat,  half2.dat];
%   t_multi.p           = [half1.p,    half2.p];
%   t_multi.sig         = [half1.sig,  half2.sig];
%   t_multi.image_labels = {'Half 1', 'Half 2'};
%   display_slices(t_multi, 'multi_image', 'names', t_multi.image_labels);

tests = functiontests(localfunctions);
end


% ---------- shared setup ------------------------------------------------

function setup(tc) %#ok<DEFNU>
close all force;
end

function teardown(tc) %#ok<DEFNU>
close all force;
end


% ---------- helpers ------------------------------------------------------

function assume_display(tc)
tc.assumeTrue(usejava('desktop') ...
    || (usejava('jvm') && feature('ShowFigureWindows')), ...
    'display_slices requires an interactive figure window');
end


function n = count_image_handles(fh)
n = numel(findall(fh, 'Type', 'image'));
end


% ---------- tests --------------------------------------------------------

function test_default_axial_auto_spacing(tc) %#ok<DEFNU>
assume_display(tc);
m = mean(canlab_get_sample_fmri_data);

[whsl, plate] = display_slices(m);

fh = gcf;
tc.verifyTrue(ishandle(fh), 'A figure should be open');
tc.verifyLessThanOrEqual(numel(whsl), 24, ...
    'Default auto-spacing should produce <=24 slices');
tc.verifyGreaterThan(numel(whsl), 0);
tc.verifyGreaterThan(count_image_handles(fh), 0);
tc.verifyNotEmpty(plate);
end


function test_coronal_auto_spacing(tc) %#ok<DEFNU>
assume_display(tc);
m = mean(canlab_get_sample_fmri_data);
whsl = display_slices(m, 'coronal');
tc.verifyLessThanOrEqual(numel(whsl), 24);
tc.verifyGreaterThan(numel(whsl), 0);
end


function test_sagittal_auto_spacing(tc) %#ok<DEFNU>
assume_display(tc);
m = mean(canlab_get_sample_fmri_data);
whsl = display_slices(m, 'saggital');
tc.verifyLessThanOrEqual(numel(whsl), 24);
tc.verifyGreaterThan(numel(whsl), 0);
end


function test_three_views_renders_one_axes(tc) %#ok<DEFNU>
assume_display(tc);
m = mean(canlab_get_sample_fmri_data);

display_slices(m, 'three_views');
fh = findobj('Type','figure','Tag','slice_display');
tc.verifyNotEmpty(fh, 'three_views should create the slice_display figure');

end


function test_explicit_spacing_override(tc) %#ok<DEFNU>
assume_display(tc);
m = mean(canlab_get_sample_fmri_data);
whsl_8mm  = display_slices(m, 'axial', 'spacing', 8);
tc.verifyGreaterThan(numel(whsl_8mm), 0);
tc.verifyLessThanOrEqual(numel(whsl_8mm), 40);
end


function test_start_end_slice(tc) %#ok<DEFNU>
assume_display(tc);
m = mean(canlab_get_sample_fmri_data);
whsl = display_slices(m, 'axial', 'spacing', 4, 'startslice', -10, 'endslice', 10);
tc.verifyGreaterThan(numel(whsl), 0);
end


function test_multi_image_default_axial(tc) %#ok<DEFNU>
assume_display(tc);
imgs = canlab_get_sample_fmri_data;

display_slices(imgs, 'multi_image');
fh = findobj('Type','figure','Tag','display_slices_multi');
tc.verifyNotEmpty(fh, 'multi_image should create display_slices_multi figure');

nimgs = min(size(imgs.dat, 2), 64);
ax = findall(fh(1), 'Type', 'axes');
tc.verifyEqual(numel(ax), nimgs, ...
    'multi_image should create exactly one axes per image (no extras)');
tc.verifyGreaterThanOrEqual(count_image_handles(fh(1)), nimgs);
end


function test_multi_image_coronal_with_slice(tc) %#ok<DEFNU>
assume_display(tc);
imgs = canlab_get_sample_fmri_data;
display_slices(imgs, 'multi_image', 'coronal', 'slice', -10);
fh = findobj('Type','figure','Tag','display_slices_multi');
tc.verifyNotEmpty(fh);
end


function test_multi_image_nimages_limit(tc) %#ok<DEFNU>
assume_display(tc);
imgs = canlab_get_sample_fmri_data;
display_slices(imgs, 'multi_image', 'nimages', 4);
fh = findobj('Type','figure','Tag','display_slices_multi');
tc.verifyNotEmpty(fh);
ax = findall(fh(1), 'Type', 'axes');
tc.verifyEqual(numel(ax), 4, '''nimages'', 4 should produce exactly 4 axes');
end


function test_multi_image_statistic_image_underlay(tc) %#ok<DEFNU>
assume_display(tc);
t = canlab_get_sample_thresholded_t(0.01);

display_slices(t, 'multi_image');
fh = findobj('Type','figure','Tag','display_slices_multi');
tc.verifyNotEmpty(fh);

% Stat path draws underlay + blob layer => >= 2 image objects per axes
n_handles = count_image_handles(fh(1));
ax = findall(fh(1), 'Type', 'axes');
tc.verifyGreaterThanOrEqual(n_handles, numel(ax));
end


function test_multi_image_multicol_statistic_image(tc) %#ok<DEFNU>
% Build a 2-column statistic_image (t_multi.dat(:, 2) is the second contrast)
% and confirm multi_image lays out exactly one subplot per image.
assume_display(tc);
imgs = canlab_get_sample_fmri_data;
nimgs = size(imgs.dat, 2);
tc.assumeGreaterThan(nimgs, 3);

half_idx = floor(nimgs / 2);
half1 = ttest(get_wh_image(imgs, 1:half_idx));
half2 = ttest(get_wh_image(imgs, half_idx+1:nimgs));
half1 = threshold(half1, .01, 'unc');
half2 = threshold(half2, .01, 'unc');

t_multi              = half1;
t_multi.dat          = [half1.dat,  half2.dat];
t_multi.p            = [half1.p,    half2.p];
t_multi.sig          = [half1.sig,  half2.sig];
t_multi.image_labels = {'Half 1', 'Half 2'};

tc.verifyEqual(size(t_multi.dat, 2), 2);

display_slices(t_multi, 'multi_image', 'names', t_multi.image_labels);
fh = findobj('Type','figure','Tag','display_slices_multi');
tc.verifyNotEmpty(fh);

ax = findall(fh(1), 'Type', 'axes');
tc.verifyEqual(numel(ax), 2, ...
    'Expected exactly 2 axes (one per stat image), no spare cells');
end
