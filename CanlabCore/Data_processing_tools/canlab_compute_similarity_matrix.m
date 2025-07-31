function [R, N] = canlab_compute_similarity_matrix(dat, varargin)
% canlab_compute_similarity_matrix: Compute pairwise similarity metrics with missing data handling
%
% - A main issue is that with brain images, zero is often treated as a missing value.
% - If elements are zero in some images and not others, correlation and other similarity metrics
% will treat these zeros as data, altering the results.
% - In a collection of images, some values (voxels) may be zero (or NaN) in some images and not others.
% - The main options are:
% 'complete cases', where all elements (voxels) must be valid for the element to be included in the similarity calculation
% 'pairwise', where only elements with valid data for each pair are included in the calculation.
%
% Pairwise is the only option that (a) uses all available data and (b)
% returns results for images a and b that do not depend on whether another
% image c is included in the set. In addition, some sets of images, e.g., neuromaps PET images, have very uneven coverage across images.
% Therefore, the pairwise missing option is generally preferred, but it is
% quite a bit slower as it requires looping through pairs (see programmers
% notes below).
%
% This function computes correlation, cosine similarity, or dot product across a set of images with pairwise missing data deletion
% It also has an option to treat zeros as data and only treat NaN values as missing.
%
% :Usage:
% ::
%     [R, N] = canlab_compute_similarity_matrix(dat, 'similarity_metric', 'cosine_similarity', 'treat_zero_as_data', true)
%
% :Inputs:
%   **dat:**
%       A [v x k] data matrix (observations x variables)
%
% :Optional Inputs:
%   **'treat_zero_as_data'** (logical):
%       Whether to include 0s as valid values. Default = false.
%
%   **'complete_cases'** (logical):
%       If true, restrict computation to rows (voxels) where all variables are valid (i.e., complete cases).
%       If false (default), compute pairwise similarity using all available data per pair.
%
%   **'similarity_metric'** (string):
%       Metric to compute. Options:
%           - 'correlation' (default)
%           - 'cosine_similarity'
%           - 'dot_product'
%           - 'dice'
%
%   **'verbose'** (logical):
%       Verbose output. Default = true.
%
%   **'doplot'** (logical):
%       Display heatmap of matrix R. Default = false.
%
% :Outputs:
%   **R:**
%       [k x k] similarity matrix
%
%   **N:**
%       [k x k] matrix of number of observations used for each pair
%
% :Author:
%   2025, Tor Wager. GPL v3.

% Programmers' Notes:
% Created by Tor Wager with GPT coding, checked/edited by Tor
% This uses a loop based approach. A non-loop approach using cov() with the 'pairwise' option was also
% explored, but this uses a different number of observations for the
% variance and cov when some elements are missing. Thus, the slower loop
% was judged to be more accurate.
% In the loop-based approach:
% 	•	For a given pair (i,j), the correlation is calculated only on rows where both x and y are non-NaN.
% 	•	The standard deviations std(x(mask)) and std(y(mask)) are computed on the same rows used for the covariance.
%
% In the cov(dat, 'partialrows') approach:
% 	•	The diagonal elements C(ii) and C(jj) represent the variance of x and y, respectively.
% 	•	But these variances are computed using all non-NaN rows in each single column, not just the shared valid rows used in the pair (i,j).

% Input parser
p = inputParser;
p.addRequired('dat', @(x) validateattributes(x, {'numeric'}, {'2d'}));
p.addParameter('treat_zero_as_data', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('similarity_metric', 'correlation', @(x) ismember(x, {'correlation', 'cosine_similarity', 'dot_product', 'dice'}));
p.addParameter('verbose', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('doplot', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('complete_cases', false, @(x) islogical(x) || isnumeric(x));  % otherwise pairwise
p.parse(dat, varargin{:});

treat_zero_as_data = logical(p.Results.treat_zero_as_data);
sim_metric = p.Results.similarity_metric;
verbose = logical(p.Results.verbose);
doplot = logical(p.Results.doplot);
complete_cases = logical(p.Results.complete_cases);

% Dice always treats zero as data
if strcmp(sim_metric, 'dice')
    treat_zero_as_data = true;
end

% Binary check for dice
if strcmp(sim_metric, 'dice')
    if ~all(ismember(dat(~isnan(dat)), [0 1]))
        error('Dice coefficient requires all input values to be 0 or 1.');
    end
end

% Replace 0s with NaN if requested
if ~treat_zero_as_data
    dat(dat == 0) = NaN;
end

% Precompute validity mask [n x k]
valid_mask = ~isnan(dat); % true if valid entry in dat

% Check whether all columns have the *same* pattern of valid values
% If so, each row in valid_mask is either all true or all false
% This will determine whether we have to loop through pairs or not
same_valid_pattern = all(valid_mask == valid_mask(:, 1), 2);

% Check if all rows are either fully valid or fully invalid
fully_consistent = all(same_valid_pattern);

wh_complete_cases = all(valid_mask, 2);

% Print result
if verbose && fully_consistent
    fprintf('All columns have consistent valid/missing value patterns. %3.0f complete cases\n', sum(wh_complete_cases));
elseif verbose
    fprintf('Columns have DIFFERENT patterns of missing or zero values. %3.0f complete cases\n', sum(wh_complete_cases));
end

if complete_cases || fully_consistent
    % No need to loop.
    dat = dat(wh_complete_cases, :);
    [R, N] = vectorized_similarity_complete_cases(dat, sim_metric);

else
        % Loop and omit missing values pairwise
    [R, N] = pairwise_corr_loop(dat, sim_metric, valid_mask);
end

[~, k] = size(dat);

% Optional plot
if doplot
    figure;
    imagesc(R); colorbar; axis square;
    title(sprintf('Similarity matrix (%s)', strrep(sim_metric, '_', '\_')));
    set(gca, 'XTick', 1:k, 'YTick', 1:k);
end

% Verbose
if verbose
    fprintf('Computed %s for %d variables. Median N = %.1f, Min N = %d\n', ...
        sim_metric, k, median(N(N>0)), min(N(N>0)));
end

end % main function





% Subfunctions
% --------------------------------------------------

function [R, N] = pairwise_corr_loop(dat, sim_metric, valid_mask)

[n, k] = size(dat);
R = nan(k);
N = zeros(k);

% Precompute pairwise combined masks [n x k x k]
% mask(:, i, j) = valid_mask(:,i) & valid_mask(:,j);
mask = false(n, k, k);
for i = 1:k
    for j = i:k
        mask(:, i, j) = valid_mask(:, i) & valid_mask(:, j);
        mask(:, j, i) = mask(:, i, j); % symmetric
    end
end

% Main loop: reuse precomputed mask
for i = 1:k
    xi = dat(:, i);
    for j = i:k
        xj = dat(:, j);
        msk = mask(:, i, j);
        Nobs = sum(msk);
        N(i, j) = Nobs;
        N(j, i) = Nobs;

        if Nobs > 1
            a = xi(msk);
            b = xj(msk);

            switch sim_metric
                case 'correlation'
                    Rval = corr(a, b);

                case 'cosine_similarity'
                    Rval = dot(a, b) / (norm(a) * norm(b));

                case 'dot_product'
                    Rval = dot(a, b);

                case 'dice'
                    Rval = 2 * sum(a & b) / (sum(a) + sum(b));

                otherwise
                    error('Unsupported similarity metric: %s', sim_metric);
            end

            R(i, j) = Rval;
            R(j, i) = Rval;
        end
    end
end


end % pairwise_corr_loop




function [R, N] = vectorized_similarity_complete_cases(dat, sim_metric)
% Efficient non-looped computation for complete cases or aligned masks

[n, k] = size(dat);
N = n * ones(k);  % All pairs have same number of observations

switch sim_metric
    case 'correlation'
        R = corr(dat);  % built-in vectorized

    case 'cosine_similarity'
        norms = sqrt(sum(dat.^2, 1));
        R = (dat' * dat) ./ (norms' * norms);  % cosine similarity

    case 'dot_product'
        R = dat' * dat;

    case 'dice'
        % For binary data only
        A = dat' * dat;             % a AND b counts
        row_sums = sum(dat, 1);     % sum of each binary map
        denom = row_sums' + row_sums;
        R = 2 * A ./ denom;

    otherwise
        error('Unsupported similarity metric: %s', sim_metric);
end

% Ensure symmetry and diagonal
R(1:k+1:end) = 1;  % optional: set diagonal to 1
end

