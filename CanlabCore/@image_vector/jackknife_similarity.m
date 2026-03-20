function [sim_values, d, low_agreement, Nvox] = jackknife_similarity(obj, varargin)
% jackknife_similarity: Jackknife spatial similarity with leave-one-image-out reference
%
% Computes similarity between each image j and the voxelwise median of the
% remaining images (N − j). This provides a robust measure of how similar
% each image is to the central tendency of the rest of the dataset.
%
% The reference image for image j is:
%
%     ref_j = median(dat(:, setdiff(1:k, j)), 2)
%
% where dat is the voxel × image matrix.
%
% Similarity is then computed between image j and ref_j using the selected
% metric.
%
% Missing data handling and options follow the conventions used in
% **canlab_compute_similarity_matrix**.
%
% :Usage:
% ::
%
%     [sim_values, d, low_agreement, Nvox] = jackknife_similarity(obj, 'similarity_metric', 'correlation')
%
% :Inputs:
%
%   **obj:**
%        An image_vector object or numeric matrix [voxels × images]
%
% :Optional Inputs:
%
%   **'similarity_metric'** (string):
%        Similarity metric used to compare each image to the median reference.
%
%        For all metrics, higher values indicate greater agreement between
%        an image and the leave-one-out group reference.
%
%        Options:
%
%        'correlation'
%            Pearson correlation coefficient between voxel values.
%            Sensitive to spatial pattern similarity after mean-centering
%            and scaling. Insensitive to overall image intensity and mean
%            shifts across voxels. Strongly influenced by relative spatial
%            patterns and voxel-wise variance across images.
%            Expected range: [-1, 1].
%
%        'cosine_similarity'
%            Cosine of the angle between voxel vectors. Equivalent to a
%            normalized dot product without mean-centering.
%            Sensitive to spatial pattern similarity and relative magnitude
%            differences across voxels. Sensitive to mean shifts and global
%            intensity changes because vectors are not centered.
%            Expected range: [-1, 1] for signed maps.
%
%        'dot_product'
%            Raw inner product between voxel vectors.
%            Sensitive to spatial pattern similarity, overall image
%            intensity, and voxel-wise magnitude. Large values may reflect
%            either similar patterns or simply large overall signal
%            magnitudes. Not normalized and therefore depends strongly on
%            voxel scaling and variance across images.
%            Expected range: unbounded (−Inf to +Inf).
%
%        'dice'
%            Dice coefficient measuring overlap between binary maps.
%            Images must be binary (0/1). Measures spatial overlap of
%            non-zero voxels and ignores voxel intensity values.
%            Insensitive to magnitude but sensitive to spatial overlap of
%            active voxels. Treats zeros as valid data.
%            Expected range: [0, 1].
%
%        'normalized_absolute_agreement'
%            Bray–Curtis similarity based on normalized L1 difference:
%
%                1 - sum(|a - b|) / sum(|a| + |b|)
%
%            Measures absolute agreement in voxel values, normalized by a metric of image intensity. 
%            Sensitive to voxel-wise magnitude differences and overall image intensity.
%            Unlike correlation, it penalizes absolute value differences
%            even when spatial patterns are proportional. Robust for sparse
%            or thresholded maps. Note: Bray-Curtis is typically used with
%            non-negative data in ecology, but can be used with signed
%            data, with this version of the denominator preferred to avoid
%            opposite-sign cancellations in the denominator.
%
%            Expected range: approximately [-1, 1] for signed maps,
%            and [0, 1] when voxel values are nonnegative.
%
%        'concordance_correlation'
%            Lin's concordance correlation coefficient (CCC). Measures both
%            correlation and agreement with the identity line a = b, and
%            therefore penalizes mean shifts, scaling differences, and
%            variance mismatches across voxels. Useful for statistical maps
%            when agreement in absolute voxel values is important.
%            Expected range: [-1, 1].
%
%   **'treat_zero_as_data'** (logical):
%        Treat zeros as valid data instead of missing values.
%        Default = false.
%
%   **'complete_cases'** (logical):
%        Restrict calculations to voxels valid for all images.
%        Default = false (pairwise deletion relative to reference).
%
%   **'verbose'** (logical):
%        Display progress and summary information.
%        Default = true.
%
%   **'doplot'** (logical):
%        Plot similarity values.
%        Default = false.
%
%   **'plot'**
%        Convenience keyword equivalent to 'doplot', true.
%        If entered, produces a diagnostic plot of jackknife similarity
%        values across images, highlighting images with low agreement.
%
%   **'verbose'** (logical)
%        Display summary information, including the group
%        agreement effect size (Cohen's d) and any images flagged as having
%        unusually low similarity to the group reference. Default = true.
%
%   **'noverbose'**
%        Convenience keyword that disables verbose output (sets 'verbose', false).
%
% :Outputs:
%
%   **sim_values:**
%        [k × 1] vector of similarity values for each image relative to
%        the median of the remaining images.
%
%   **d:**
%        Cohen's d effect size for the group. Higher values indicate good group
%        agreement. The expected effect size under random noise is 0.
%        Note this is not a true Cohen's d for underlying metrics that
%        are bounded, and is in these cases a "pseudo-d".
%
%   **low_agreement:**
%        Logical vector of images whose jackknife similarity is >3 MAD below the median,
%        indicating unusually low similarity to the group reference based on other images.
%        If images are expected to be similar, these images can be
%        considered outliers.
%
%   **Nvox:**
%        [k × 1] number of voxels used in each similarity calculation.
%
% :Recommendations on using metrics to evaluate participant-level
% contrast images subjected to group analysis:
%
% For a single metric, standardized absolute deviation is the most
% appealing as a single metric for agreement across independent
% contrast images. It is sensitive to pattern noise and random
% shifts in location or scale of individual contrast images.
% Sensitivity to all three metrics is preserved in the group effect
% sizes (Cohen''s d). This metric is best interpreted as a
% normalized L1 deviation from the out-of-sample group reference,
% analogous to the average absolute residual in a voxelwise group
% model. This metric is also most similar to the error in a
% voxel-wise group analysis, averaged across voxels.
%
% Cosine similarity is our second choice for agreement across
% independent contrast images (i.e., participants). Like
% correlation, it is sensitive to pattern noise in images, and also
% sensitive to random shifts in the image location (global
% activation or deactivation) in a way that correlation is not.
% Note, however, that it is sensitive to changes in the location of
% the entire image, i.e., the level of global activation or
% deactivation across an entire set of contrast images, and global
% activation in contrast images can increase cosine similarity.
% This property introduces interpretational ambiguity, as higher
% values may be a function of whole-dataset global values rather
% than agreement. Like correlation, the absolute magnitude is
% interpretable, and unlike correlation, group effect sizes
% (Cohen''s d) within a set of contrast images remain sensitive to
% both pattern noise and random shifts in image location. Thus, the
% overall d value — which reflects both the mean level of agreement
% and the spread across individual images — can be interpreted,
% with higher values reflecting greater inter-image agreement.
%
% Correlation is also sensible as a single metric for agreement
% across independent contrast images. It is sensitive to pattern
% noise, but not to global shifts in location or scale or random
% shifts in location or scale of individual contrast images. The
% former two properties are a strength, whereas the latter two are
% a drawback, as agreement should decrease with these types of
% random shifts. In addition, though the magnitude of correlation
% is interpretable, group effect sizes are surprisingly insensitive
% even to pattern noise across the entire dataset. Because
% correlation removes mean and scale, global pattern noise can
% reduce similarity uniformly across all images. This reduces
% within-group variance relative to the mean similarity,
% artificially inflating effect size (d). Thus, correlation group
% effect size should be used with caution.
%
% As a suite of metrics, correlation, mean shift z and scale shift
% z form a coherent set. Each is preferentially sensitive to a
% different type of inter-participant noise. Correlation is
% sensitive to pattern noise, mean shift z to inter-image variation
% in location (global mean), and scale shift z to to inter-image
% variation in scale (multiplicative noise factors). They thus
% provide an orthogonal decomposition of agreement into pattern
% similarity (correlation), location shifts (mean_shift_z), and
% scale shifts (scale_shift_z). Both mean shift z and scale shift z
% are sensitive to these noise sources in both mean raw values and
% group effect sizes (d).
%
% Other measures have both some desirable properties and some
% drawbacks. Dot product and normalized absolute agreement
% (Bray-Curtis) are sensitive to global whole-dataset location
% shifts and have no superior properties to other metrics. Lin's
% Concordance Correlation performs most similiarly to cosine
% similarity, sensitive to random image-level location shifts and
% pattern noise, but like correlation has the drawback of
% insensitivity to pattern noise at the whole-dataset effect size
% level.
%
% :Examples:
% ::
%
%    % Load CANlab sample images
%    obj = load_image_set('emotionreg');
%
%    % Compute jackknife correlation similarity
%    [sim_values, Nvox] = jackknife_similarity(obj);
%
%    % Use cosine similarity instead
%    sim_cos = jackknife_similarity(obj, 'similarity_metric', 'cosine_similarity');
%
%    % Plot results
%    jackknife_similarity(obj, 'doplot', true);
%
%    % Treat zeros as data
%    sim_values = jackknife_similarity(obj, 'treat_zero_as_data', true);
%
%
% :Author:
%   2026, Tor Wager. GPL v3.
%
<<<<<<< Updated upstream
=======
% See also:
% jackknife_similarity_unit_test.m

>>>>>>> Stashed changes
% Programmers' Notes:
% This function follows conventions used in canlab_compute_similarity_matrix.
% The jackknife reference is recomputed for each image to avoid circularity.
%

% -------------------------------------------------------------------------
% Input parsing
% -------------------------------------------------------------------------

% Parse special command keywords and remove them before inputParser
doplot = false;  
plot_idx = strcmpi(varargin, 'plot');
if any(plot_idx)               % Override: omit 'doplot' key/value pair 
    doplot = true;
    varargin(plot_idx) = [];   % remove so inputParser doesn't see it
end

verbose = true;  
verbose_idx = strcmpi(varargin, 'noverbose');
if any(verbose_idx)
    verbose = false;
    varargin(verbose_idx) = [];   % remove so inputParser doesn't see it
end

% Use inputParser to parse key/value pairs
% First add obligatory/non-conditional keywords
p = inputParser;
p.addRequired('obj');
p.addParameter('similarity_metric','correlation',...
<<<<<<< Updated upstream
    @(x) ismember(x,{'correlation','cosine_similarity','dot_product','dice','absolute_agreement'}));
=======
    @(x) ismember(x, {'correlation','cosine_similarity','dot_product', ...
                      'dice','normalized_absolute_agreement','concordance_correlation', ...
                      'standardized_abs_deviation','mean_shift_z','scale_shift_z'}));
>>>>>>> Stashed changes
p.addParameter('treat_zero_as_data', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('complete_cases', false, @(x) islogical(x) || isnumeric(x));

% Special key/value pairs that we have potentially set with optional keywords
p.addParameter('doplot', doplot, @(x) islogical(x) || isnumeric(x));
p.addParameter('verbose', verbose, @(x) islogical(x) || isnumeric(x));

% process inputs and deal out to variables in workspace
p.parse(obj,varargin{:});

sim_metric = p.Results.similarity_metric;
treat_zero_as_data = logical(p.Results.treat_zero_as_data);
complete_cases = logical(p.Results.complete_cases);
verbose = logical(p.Results.verbose);
doplot = logical(p.Results.doplot);

% -------------------------------------------------------------------------
% Initial checks and setup
% -------------------------------------------------------------------------

if isa(obj, 'image_vector')
    dat = obj.dat;
else
    dat = obj;
end

<<<<<<< Updated upstream
=======
% -------------------------------------------------------------------------
% Preprocess
% -------------------------------------------------------------------------

% Data integrity check: ensure no empty / zero-norm images
img_L1 = sum(abs(dat), 1, 'omitnan');

if any(~isfinite(img_L1)) || any(img_L1 == 0)
    error(['jackknife_similarity: invalid images detected. ' ...
        'Some images have zero L1 norm or contain only NaN/zero values.']);
end

>>>>>>> Stashed changes
% Dice requires binary data
if strcmp(sim_metric,'dice')

    treat_zero_as_data = true;

    dat_nonan = dat(~isnan(dat));

    if any(dat_nonan ~= 0 & dat_nonan ~= 1)
        error('Dice coefficient requires binary maps coded as 0/1.');
    end

end

if ~treat_zero_as_data
    dat(dat==0) = NaN;
end

[nvox, k] = size(dat);

sim_values = nan(k,1);
Nvox = zeros(k,1);

% -------------------------------------------------------------------------
<<<<<<< Updated upstream
=======
% Stability constant and summaries
% -------------------------------------------------------------------------

% dat_valid = dat(~isnan(dat));
% global_abs_scale = median(abs(dat_valid));
% if isempty(global_abs_scale) || global_abs_scale == 0
%     global_abs_scale = 1;
% end
% epsilon = eps; % 1e-6 * global_abs_scale;

% The role for epsilon:
% prevent division by zero or numerical explosions when across-image spread is
% nearly zero, without materially changing ordinary cases. A tiny fraction of
% the global absolute scale is more sensible than an arbitrary constant like
% eps or 1e-8, because it stays in the same measurement units as the images.

if strcmpi(sim_metric, 'scale_shift_z') || strcmpi(sim_metric, 'mean_shift_z') || strcmpi(sim_metric, 'standardized_abs_deviation')

    % robust within-image scale after centering
    img_scale = mad(dat);
    global_scale = median(img_scale);

    img_mean = mean(dat,1,'omitnan')';
    % 
    % % Global normalization scales (do NOT depend on across-image variability)
    % 
    % img_centered = dat - img_mean';
    % img_scale = median(abs(img_centered),1,'omitnan')';
    % 
    % global_scale = mad(img_scale);   % stable reference scale
    % 
    % % % also useful for mean shift normalization:
    % % global_mean_scale = mean(abs(img_mean), 'omitnan');

end


% -------------------------------------------------------------------------
>>>>>>> Stashed changes
% Main loop
% -------------------------------------------------------------------------

idx = true(1, k);

for j = 1:k

    % Faster to use a logical vector we flip bits on than setdiff
    % others = setdiff(1:k, j);
    others = idx;
    others(j) = false;

    ref = median(dat(:, others), 2, 'omitnan');
    x = dat(:, j);

    if treat_zero_as_data
        % zeros are treated as data

        valid = ~isnan(x) & ~isnan(ref);

        if complete_cases
            valid = valid & all(~isnan(dat(:, others)), 2);
        end

    else
        % zeros are treated as missing

        valid = ~isnan(x) & x ~= 0 & ~isnan(ref) & ref ~= 0;

        if complete_cases
            valid = valid & all(~isnan(dat(:, others)), 2) & all(dat(:, others) ~= 0, 2);
        end

    end

    a = x(valid);
    b = ref(valid);

    Nvox(j) = numel(a);

    if Nvox(j) < 2
        sim_values(j) = NaN;
        continue
    end

<<<<<<< Updated upstream
=======
    % if strcmp(sim_metric,'standardized_abs_deviation')
    %     ref_scale = row_mad_omitnan(dat(:, others));
    %     s = ref_scale(valid);
    % else
    %     ref_scale = [];
    % end


>>>>>>> Stashed changes
    switch sim_metric

        case 'correlation'

            sim_values(j) = corr(a,b);


        case 'cosine_similarity'
            % sim_values(j) = dot(a,b) / (norm(a) * norm(b));

            den = norm(a) * norm(b);
            if den == 0
                sim_values(j) = NaN;  % return NaN if either norm = 0
            else
                sim_values(j) = dot(a,b) / den;
            end

        case 'dot_product'
            sim_values(j) = dot(a,b);

        case 'dice'

            % a = a ~= 0;  % enforce exactly 0/1 values (though checked above)
            % b = b ~= 0;
            sim_values(j) = 2 * sum(a & b) / (sum(a) + sum(b));

        case 'normalized_absolute_agreement'
            % sim_values(j) = 1 - sum(abs(a - b)) ./ sum(abs(a + b));
            % sim_values(j) = 1 - sum(abs(a - b)) / sum(abs(a) + abs(b)); % different than above for signed data

            den = sum(abs(a) + abs(b));
            if den == 0
                sim_values(j) = NaN;
            else
                sim_values(j) = 1 - sum(abs(a - b)) ./ den;
            end

        case 'concordance_correlation'

            ma = mean(a);
            mb = mean(b);

            va = var(a, 1);   % population variance
            vb = var(b, 1);
            cab = mean((a - ma) .* (b - mb));

            denom = va + vb + (ma - mb)^2;

            if denom == 0
                sim_values(j) = NaN;
            else
                sim_values(j) = 2 * cab / denom;
            end

<<<<<<< Updated upstream
=======
        case 'standardized_abs_deviation'

            % Mean absolute voxelwise deviation, standardized by leave-one-out
            % voxelwise MAD, then transformed so higher = better agreement.
            shift = mean(abs(a - b) ./ global_scale);
            sim_values(j) = 1 ./ (1 + shift);

        case 'mean_shift_z'

            others_mu = img_mean(others);
            mu_center = mean(others_mu, 'omitnan');
            % mu_spread = scalar_std_omitnan(others_mu);

            shift = abs(img_mean(j) - mu_center) ./ global_scale; %(mu_spread + epsilon);
            
            sim_values(j) = 1 ./ (1 + shift);

        case 'scale_shift_z'

            others_scale = img_scale(others);
            scale_center = mean(others_scale, 'omitnan');
            % scale_spread = scalar_std_omitnan(others_scale);

            shift = abs(img_scale(j) - scale_center) ./ global_scale; % (scale_spread + epsilon);
            sim_values(j) = 1 ./ (1 + shift);

>>>>>>> Stashed changes
        otherwise
            error('Unsupported similarity metric.');
    end

end

% d: Cohen's d effect size for the group. Values > 0 indicate good group
% agreement. The expected effect size under random noise is 0.
if std(sim_values) == 0
    d = Inf;
else
    d = mean(sim_values) ./ std(sim_values);
end

% low_agreement: Logical vector of images whose jackknife similarity is >3 MAD below the median,
% indicating unusually low similarity to the group reference based on other images.
% low_agreement = (sim_values - (median(sim_values) - (3 * mad(sim_values)))) < 0;
low_agreement = sim_values < median(sim_values) - 3 * mad(sim_values);

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------

if doplot

    create_figure('jackknife_similarity');

    plot(sim_values, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', [0.7 0.7 0.7]); 
    hold on

    % highlight low-agreement images
    if any(low_agreement)
        plot(find(low_agreement), sim_values(low_agreement), 'ro', ...
            'MarkerSize', 10, 'LineWidth', 2);
    end

    % reference threshold line
    thr = median(sim_values) - 3 * mad(sim_values);
    yline(thr, '--', '3 MAD threshold', 'LineWidth', 1.5);

    xlabel('Image');
    ylabel(['Similarity (' strrep(sim_metric,'_',' ') ')']);

    titlestr = sprintf('Jackknife %s (d = %.2g)', strrep(sim_metric, '_', ' '), d); 
    title(titlestr);

    box on

end

% -------------------------------------------------------------------------
% Verbose
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Verbose output
% -------------------------------------------------------------------------

if verbose

    fprintf('\nJackknife similarity computed for %d images\n', k);
    fprintf('Similarity metric: %s\n', sim_metric);
    fprintf('Median number of voxels used: %.0f\n', median(Nvox));
    fprintf('Group agreement effect size (Cohen''s d): %.2g\n', d);

    if any(low_agreement)

        idx = find(low_agreement);

        fprintf('Images with low agreement (>3 MAD below median):\n');
        fprintf('  %s\n', num2str(idx'));

    else

        fprintf('No low-agreement images detected.\n');

    end

end

<<<<<<< Updated upstream
end % jackknife_similarity
=======
end % jackknife_similarity


% -------------------------------------------------------------------------
% Helpers
% -------------------------------------------------------------------------

function m = row_mad_omitnan(X)
medx = median(X,2,'omitnan');
m = median(abs(X - medx),2,'omitnan');
end

% function s = scalar_std_omitnan(x)
% x = x(~isnan(x));
% if numel(x)<2
%     s = NaN;
% else
%     s = std(x);
% end
% end
>>>>>>> Stashed changes
