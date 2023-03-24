function [mad_similarity cos_similarity r_similarity] = outliers_xval(obj)
% Similarity metrics comparing each image to the mean of others in the set using repeated k-fold cross-validation
%
% :Usage:
% ::
%
%     [mad_similarity cos_similarity r_similarity] = outliers_xval(obj)
%
% :Background:
% ::
% Assessing image outliers is crucial for detecting artifacts and other
% violations of assumptions underlying image analysis. 
% 
% Many methods are available, but they often use measures of variance defined relative to a sample,
% which creates problems. For example, Mahalanobis distance uses squared
% multivariate distances among images. This is useful, but detects only
% relative outliers. If there are many corrupted images, all images will
% look normal by comparison to others.  Inter-image correlations are also
% commonly used. These are also useful, but they are also sensitive to the
% overall distribution of artifacts across the set. If there are many
% corrupted images, correlations will tend to be low overall, and the bad 
% images will be harder to detect. In cases where comparing individual
% images to a gold standard is desired, but there is no external standard
% image (perhaps because images look somewhat different in each sample),
% comparing images to a mean image derived from others in the set may be useful. 
% In this case, however, it's desirable for the mean to be derived
% independent of the test image itself. This can be accomplished with
% cross-validation, bootstrapping, or jackknife analyses.
%
% Here, we compare each of a set of images to a mean derived from other
% images in the set using k-fold cross-validation. To average over errors
% related to the particular fold splits chosen, we average over a number of
% repeated cross-valiations.
%
% We return two error metrics, which have different desirable properties
% depending on the type of image and use case.
% - Pearson's correlation (r): correlation is insensitive to mean shift and scale
%   it is most useful when only the pattern in the image is meaningful
%   and/or images will be rescaled in analysis, and the zero-point and
%   scale are not important.
%
% - cosine similarity is sensitive to mean shift, not scale
%   This is useful for images with a meaningful zero-point that should match across images,
%   but where the scale is irrelevant (or will be removed, e.g., via
%   normalization)
%
% - Absolute agreement is most relevant when mean and scale are both relevant
%   This is the case with many images, e.g., those subjected to group statistical analysis of intensity
%   values. e.g., when testing task - control contrast images across a group,
%   We are interested in whether the group mean is 0. The zero-point is meaningful and 
%   it should agree if the individual test images are replicates of the same effect.
%   The scale is also meaningful, as deviations in the values determine the variance estimate 
%   of the test statistic. Any deviation is treated as error, even if it results from 
%   mean shift or scale varation in the image as a whole.  
%
%
% :Inputs:
%
%   **obj:**
%        an image_vector object (e.g., fmri_data) with >=5 images
%
% :Outputs:
%
%   **mad_similarity:**
%        [nimages x 1] vector of absolute agreement between each image and the mean, 
%        1/(1+MAD), averaged across repeated k-fold
%
%   **cos_similarity:**
%        [nimages x 1] vector of cosine similarity between each image and the mean, 
%        averaged across repeated k-fold
%
%   **r_similarity:**
%        [nimages x 1] vector of correlation between each image and the mean, 
%        averaged across repeated k-fold
%
% ----------------------------------------------------------------------
% Examples:
% ----------------------------------------------------------------------
% test_images = load_image_set('emotionreg');
%
% [mad_similarity cos_similarity r_similarity ] = outliers_xval(test_images);
% figure; plot(mad_similarity)
% hold on;
% plot(cos_similarity)
% plot(r_similarity)
% legend({'Abs agreement' 'Cos sim' 'r'})
%
% % Compare with outliers based on Mahalanobis and other metrics
% figure;  
% [est_outliers_uncorr, est_outliers_corr, outlier_tables] = outliers(test_images, 'notimeseries');
%
% Another, larger sample dataset
% test_images = load_image_set('kragel18_alldata');
% [mad_similarity cos_similarity r_similarity ] = outliers_xval(test_images);

% ..
%     Author and copyright information:
%
%     Copyright (C) 2023 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..

% ----------------------------------------------------------------------
% Default parameters
% ----------------------------------------------------------------------

nfolds = 5;
nrepeats = 100;  % number of cross-validation repeats

% ----------------------------------------------------------------------
% Set up variables and functions
% ----------------------------------------------------------------------

% Analyze gray matter, where we expect meaningful variation
%obj_orig = obj;
obj = apply_mask(obj, fmri_data(which('gray_matter_mask.nii'), 'noverbose'));

nimages = size(obj.dat, 2);

if nimages < 5, error('Object must contain at least 5 images'), end

Yid = ones(nimages, 1);

[r_test, cos_sim_test, mad_test] = deal(NaN .* zeros(nimages, nrepeats));

cvpart = cvpartition(Yid, 'k', nfolds);
% [S.trIdx, S.teIdx] = xval_stratified_holdout_leave_whole_subject_out(S.Y, S.id, 'doverbose', doverbose, 'doplot', doplot, 'nfolds', nfolds);

cos_sim =  @(x, y) x' * y ./ (norm(x) * norm(y));

% ----------------------------------------------------------------------
% Run
% ----------------------------------------------------------------------

for p = 1:nrepeats

    % Draw a new k-fold cross-validation partition
    cvpart = cvpart.repartition;

    for i = 1:nfolds

        train_set = get_wh_image(obj, cvpart.training(i));
        test_set = get_wh_image(obj, cvpart.test(i));

        m = mean(train_set);

        % correlation is insensitive to mean shift and scale
        % cosine sim is sensitive to mean shift, not scale
        % for images with a meaningful zero-point that should match across images,
        % cosine sim is more meaningful.
        % if mean and scale are both relevant, as with contrast images, absolute
        % agreement may be the best metric, e.g. MAD = median absolute deviation

        r_test(cvpart.test(i), p) = corr(m.dat, test_set.dat)';

        cos_sim_test(cvpart.test(i), p) = cos_sim(m.dat, test_set.dat)';

        mad_test(cvpart.test(i), p) = median(abs(m.dat - test_set.dat))';

    end % k-fold

end % cv repeats

% ----------------------------------------------------------------------
% Calculate final similarity measures
% ----------------------------------------------------------------------

r_similarity = mean(r_test')';
cos_similarity = mean(cos_sim_test')';

mad_dissim = mean(mad_test')'; % this is actually still deviation here
mad_similarity = 1./(1+mad_dissim); % convert to similarity, scale from [0 1] 
% +1 scales so that zero error (perfect agreement) will have a value of 1.

end % function

