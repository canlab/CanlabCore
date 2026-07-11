function tests = canlab_test_resampling_pattern_expression
%CANLAB_TEST_RESAMPLING_PATTERN_EXPRESSION Pattern expression is stable across resampling.
%
% apply_mask(..., 'pattern_expression') resamples the weight map into the
% data space when the two differ. This test checks that the resulting
% pattern-expression scores barely change regardless of which interpolation
% method does that resampling (default vs nearest vs spline) - i.e. the
% dot-product readout is not an artifact of the interpolation choice.
%
% Uses the SIIPS signature against the emotionreg sample. Skipped if the
% SIIPS image set is not available on the path (it ships with
% Neuroimaging_Pattern_Masks).
%
% Converted from the old standalone visual script
% Unit_tests/old_to_integrate/resampling_pattern_expression_unit_test1.m,
% which only plotted the variants and asserted nothing.

tests = functiontests(localfunctions);
end


function test_pattern_expression_stable_across_interp(tc)   %#ok<*DEFNU>
obj = canlab_get_sample_fmri_data();

try
    siips = load_image_set('siips', 'noverbose');
catch ME
    tc.assumeFail(['SIIPS signature not available on this runner: ' ME.message]);
end
tc.assumeNotEmpty(siips.dat, 'SIIPS loaded but empty');

n = size(obj.dat, 2);

% Default path lets apply_mask resample internally; the other two pre-resample
% the weight map with an explicit interpolation method.
pe_default = apply_mask(obj, siips, 'pattern_expression');
pe_nearest = apply_mask(obj, resample_space(siips, obj, 'nearest'), 'pattern_expression');
pe_spline  = apply_mask(obj, resample_space(siips, obj, 'spline'),  'pattern_expression');

for v = {pe_default, pe_nearest, pe_spline}
    tc.verifyEqual(numel(v{1}), n);
    tc.verifyTrue(all(isfinite(v{1})), 'pattern expression returned non-finite values');
end

% Empirically these correlate at >= 0.9999; 0.99 is a safe regression floor.
tc.verifyGreaterThan(corr(pe_default, pe_nearest), 0.99);
tc.verifyGreaterThan(corr(pe_default, pe_spline),  0.99);
end
