function t = canlab_get_sample_thresholded_t(p_thresh)
%CANLAB_GET_SAMPLE_THRESHOLDED_T Build a thresholded statistic_image for tests.
%
%   t = canlab_get_sample_thresholded_t()
%   t = canlab_get_sample_thresholded_t(p_thresh)   % default 0.05
%
%   Loads the emotionreg sample, runs a one-sample voxelwise t-test, and
%   thresholds the resulting statistic_image at p < p_thresh (uncorrected).
%   Used by display / table / extraction tests as a shared input fixture
%   so each one doesn't recompute the t-test independently.

if nargin < 1, p_thresh = 0.05; end

obj = canlab_get_sample_fmri_data();
t = ttest(obj);
t = threshold(t, p_thresh, 'unc');
end
