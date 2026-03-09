function resampling_pattern_expression_unit_test1

obj = load_image_set('emotionreg'); % test dataset
siips = load_image_set('siips');    % sampled to higher-res mask on loading
siips_resamp = resample_space(siips, obj, 'spline');
siips_orig = fmri_data(siips.image_names); % original pattern

val = [];
val(:, 1) = apply_mask(obj,  siips_orig, 'pattern_expression');
val(:, 2) = apply_mask(obj,  resample_space(siips_orig, obj), 'pattern_expression');
val(:, 3) = apply_mask(obj,  resample_space(siips_orig, obj,  'nearest'), 'pattern_expression');
val(:, 4) = apply_mask(obj,  resample_space(siips_orig, obj,  'spline'), 'pattern_expression'); 

val(:, 5) = apply_mask(obj,  siips, 'pattern_expression');
val(:, 6) = apply_mask(obj,  resample_space(siips, obj), 'pattern_expression');
val(:, 7) = apply_mask(obj,  resample_space(siips, obj,  'nearest'), 'pattern_expression');
val(:, 8) = apply_mask(obj,  resample_space(siips, obj,  'spline'), 'pattern_expression');

val(:, 9) = apply_mask(obj,  siips, 'pattern_expression');
val(:, 10) = apply_mask(obj,  resample_space(siips, obj), 'pattern_expression');
val(:, 11) = apply_mask(obj,  resample_space(siips, obj,  'nearest'), 'pattern_expression');
val(:, 12) = apply_mask(obj,  resample_space(siips, obj,  'spline'), 'pattern_expression');

val(:, 13) = apply_mask(obj,  siips_resamp, 'pattern_expression');
sdat = apply_siips(obj);
val(:, 14) = sdat{1};

m = mean(val);
s = std(val);
c = corr(val);

create_figure('siips comparison', 1, 3); 
errorbar(m, s); title('means and standard errors')
subplot(1, 3, 2); plot(m); title('means (zooming in on differences')
subplot(1, 3, 3); imagesc(c); colorbar
title('Correlations among variants')

set(gca, 'YDir', 'reverse'); axis tight

end
