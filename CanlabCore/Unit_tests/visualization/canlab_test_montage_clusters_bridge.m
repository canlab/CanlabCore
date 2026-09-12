function tests = canlab_test_montage_clusters_bridge
%CANLAB_TEST_MONTAGE_CLUSTERS_BRIDGE Tests for canlab_montage_clusters and
%   the fmridisplay-based replacements of the legacy montage_clusters
%   branches in region.montage, image_vector.montage and fmri_data.plot.
tests = functiontests(localfunctions);
end


function setup(tc)
tc.assumeTrue(usejava('jvm'), 'Requires a MATLAB session with graphics (JVM).');
close all
end


function teardown(~)
close all
end


function r = get_regions()
persistent rcache
if isempty(rcache)
    t = canlab_get_sample_thresholded_t();
    evalc('rcache = region(t);');
end
r = rcache;
end


% -------------------------------------------------------------------------
% canlab_montage_clusters
% -------------------------------------------------------------------------

function test_region_input_returns_fmridisplay(tc)
r = get_regions();
o2 = canlab_montage_clusters(r, 'colors', {[1 .3 0]}, 'noverbose', 'figname', 'bridge test');
tc.verifyClass(o2, 'fmridisplay');
tc.verifyEqual(numel(o2.activation_maps), 1);
tc.verifyEqual(numel(o2.activation_maps{1}.blobhandles) > 0, true);
fh = findall(0, 'Type', 'figure', 'Name', 'bridge test');
tc.verifyEqual(numel(fh), 1, 'figname should be applied to the montage figure');
end


function test_struct_input_sagittal(tc)
% Legacy clusters structure on sagittal slices (montage_clusters_medial replacement)
r = get_regions();
cl = region2struct(r);
o2 = canlab_montage_clusters(cl, 'colors', {'r'}, 'montagetype', 'sagittal', 'noverbose');
tc.verifyClass(o2, 'fmridisplay');
tc.verifyEqual(numel(o2.activation_maps), 1);
tc.verifyEqual(o2.montage{1}.orientation, 'sagittal');
end


function test_multiple_layers_and_empty_layer(tc)
r = get_regions();
tc.assumeTrue(numel(r) >= 2, 'Need at least two regions');
empty_struct = region2struct(r); empty_struct = empty_struct([]);   % 0-element struct array
o2 = canlab_montage_clusters({r(1), r(2:end), empty_struct, []}, 'colors', {'g' 'b'}, 'noverbose');
tc.verifyEqual(numel(o2.activation_maps), 2, 'Empty layers must be skipped');
end


function test_colormap_flag_and_existing_o2(tc)
r = get_regions();
o2 = canlab_results_fmridisplay(region(), 'compact2', 'noblobs', 'nooutline', 'noverbose');
fh_before = findall(0, 'Type', 'figure');
o2 = canlab_montage_clusters(r, 'colormap', 'o2', o2, 'noverbose');
tc.verifyEqual(numel(o2.activation_maps), 1);
fh_after = findall(0, 'Type', 'figure');
tc.verifyEqual(numel(setdiff(fh_after, fh_before)), 0, 'Passing o2 must not open a new montage figure');
end


function test_bad_layer_errors(tc)
tc.verifyError(@() canlab_montage_clusters({42}, 'noverbose'), 'canlab_montage_clusters:BadInput');
end


% -------------------------------------------------------------------------
% Migrated legacy branches
% -------------------------------------------------------------------------

function test_region_montage_old_option_uses_fmridisplay(tc)
r = get_regions();
tc.verifyWarning(@() montage(r, 'old', 'noverbose'), 'region:montage:OldDeprecated');
end


function test_image_vector_montage_scnmontage(tc)
t = canlab_get_sample_thresholded_t();
obj = fmri_data(); obj.volInfo = t.volInfo; obj.removed_voxels = t.removed_voxels;
obj.dat = [t.dat .* double(t.sig), -t.dat .* double(t.sig)];
fh = montage(obj, 'scnmontage');
tc.verifyEqual(numel(fh), 2);
tc.verifyTrue(all(isgraphics(fh, 'figure')));
tc.verifyEqual(get(fh(2), 'Tag'), 'Montage   2');
end


function test_fmri_data_plot_means_for_unique_Y(tc)
obj = canlab_get_sample_fmri_data();
obj.Y = [ones(15, 1); 2 * ones(15, 1)];
evalc('plot(obj, ''means_for_unique_Y'');');
fh = findall(0, 'Type', 'figure', 'Name', 'Montage_mean_across_conditions');
tc.verifyEqual(numel(fh), 1);
fh = findall(0, 'Type', 'figure', 'Name', 'Montage_coeff_of_var_across_conditions');
tc.verifyEqual(numel(fh), 1);
end


function test_montage_clusters_deprecation_warning(tc)
r = get_regions();
cl = region2struct(r(1));
clear montage_clusters   % reset the once-per-session flag
tc.verifyWarning(@() montage_clusters([], cl, {'r'}), 'montage_clusters:Deprecated');
end
