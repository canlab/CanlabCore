function similarity_output = canlab_pattern_similarity(dat, pattern_weights, varargin)
% Calculate similarity between each column in a data matrix dat and a vector of pattern weights
%
% - Similarity options: dot product, cosine similarity, and correlation
% - Columns are often images, e.g., from fmri_data.dat for fmri_data objects
% - weights are pattern weights from one or more 'signature' patterns or binary masks (each pattern is a column)
%
% - Assumes dat and pattern_weights matrices have equal rows (voxels) and include valid voxels
% - Removes empty column data column-wise (image-wise) in case some images have uneven voxel coverage
%
% - Used in apply_mask, extract_roi_averages, image_similarity_plot***
%
% similarity_output = canlab_pattern_similarity(dat, pattern_weights)
%
% :Inputs:
%   **dat**  voxels x images matrix of data to compare to patterns
%
%   **pattern_weights** voxels x images matrix of patterns to compare to data
%
% :Optional Inputs:
%
%   **dot_product**
%       [Default] Use dot product
%
%   **cosine_similarity**
%       Use cosine similarity metric for pattern expression instead of dot product
%
%   **correlation, corr**
%       Use correlation metric for pattern expression instead of dot product
%
%   **binary_overlap**
%       Use percent overlap of binary masks instead of dot product. needs
%       binary masks and pattern weights.
%
%   **posterior_overlap**
%       Using percent overlap of binary masks, this option calculates the
%       posterior probability of observing non-zero pattern weights (binary)
%       given binary masks
%
%   **ignore_missing**
%       Suppress printing of warnings when thre are missing values
%       (zeros/NaN) in data images. This function always removes voxels
%       from analysis with missing data values. This does not affect the
%       dot product metric, but does affect cosine similarity and
%       correlation.
%
%   **no_warnings**
%       Suppress output on missing values.
%
%
%   **treat_zero_as_data**
%       The default behavior is to calculate image similarity on complete cases pairwise. 
%       Complete cases means a voxel has a valid value (non-zero, non-NaN) for both images being correlated. 
%       Otherwise, values that are zero in one image and non-zero and another will be 
%       part of the correlation and will have a strong influence on the correlation value.
%       Zero is treated as a missing value by convention (as per SPM
%       standards).
%
%       In some situations, e.g., when correlating binary masks (1/0 for both) or correlating 
%       a binary mask with another image, the zero value within data.obj could be
%       meaningful. Here data.obj is a thresholded map (0 means value under-threshold 
%       rather than being a missing value) or binary map. In this case, you
%       can use 'treat_zero_as_data', 1 to treat the zeros as data values.
%
%   **exclude_zero_mask_values**
%
% :Outputs:
%   **similarity_output**
%       Matrix of similarity measures.
%
%
% :Examples:
% ::
% dat = rand(100, 5);  % 100 voxels, 5 images
% weights = rand(100, 2); % 100 voxels, 2 patterns
% similarity_output = canlab_pattern_similarity(dat, weights, 'cosine_similarity');
%
% :See also:
%
% apply_mask
% image_similarity_plot
% extract_roi_averages, to get individual region averages / local pattern expression
% apply_nps, which does whole-pattern and local regional expression
%
% ..
%    Notes:
%
% Dealing with missing data and partial coverage
% ---------------------------------------------------
% Pattern pattern_weights can have zeros, which may be valid values in voxels,
% i.e., with binary masks
% By default, data images with values of zero or NaN are considered out-of-mask, as they
% are not valid values. That is, 0 is typically treated as a missing data value in images
% This can be changed by setting ['treat_zero_as_data', true].  
% This function does not exlude 0 voxels from the pattern_mask input, where it is assumed 
% that 0 is a meaninful value (for a neuromarker or binary mask). 
% Thus, this function treats values of 0 in DATA images, not pattern masks,
% as missing values, and excludes these voxels from both image and
% mask when calculating similarity.
% Thus, there is an asymmetry between pattern mask and image data
% in considering which voxels to use.
% Otherwise, all dot product, correlation, and other similarity metrics are standard.
%
% When comparing two sets of binary images (e.g. k-means clusters) and the
% cluster of the input solution does not overlap with any of the target
% patterns/clusters, cosine similarity is attempting division by 0. Instead
% of returning NaN, we set those cosine values to 0.
%
% Effects of zeros/missing values on similarity metrics
% ---------------------------------------------------
% Dot product is affected by voxel size, coverage in image and mask, scale
% in general.
% Cosine similarity is affected by coverage in image
% If we remove voxels from analysis that are not in data image (NaN or zero) first,
% then cosine similarity is unaffected by coverage, in the sense that the
% upper bound is 1.
% Correlation is affected by coverage in image and mean level in image. If
% there are zeros in image, mean-centering will be affected.
% If we remove NaN/zero values from analysis on an image-by-image basis, we
% are computing partial similarity for the areas of the image that exist.
% This is good in terms of preserving scale and analyzing the data we have,
% but changes the measurement properties of the pattern measures, because
% we're looking only at a subset of the model.
%
% We also need to calculate bad values on an image-by-image basis, not
% relying on remove_empty to exclude voxels with ineligible values
% across the entire set of images.
%
%  Created by Tor Wager - 3/7/17

% Programmer's Notes:
%
% 2017/09/07 Stephan Geuter
%   - added option for percent overlap of binary masks (see also
%   image_similarity_plot.m and riverplot.m
%   - Changed metric selection to string format.
%
% 2022/03/31 Ke Bo
%   - added option for treating zero value in the map as real value rather
%   than missing data
%

% ---------------------------------
% Defaults and optional inputs
% ---------------------------------

% docorr = false;     % run correlation instead of dot-product for pattern expression
sim_metric = 'dotproduct'; % Default: Correlation. SG. docosine = false;   % run cosine sim instead
doignoremissing = false; % ignore warnings for missing voxels
doprintwarnings = true;  % print warnings regarding missing voxels, etc.
treat_zero_as_data = false; % Treat zero value as missing data.
exclude_zero_mask_values = false;

if any(strcmp(varargin, 'ignore_missing'))
    doignoremissing = true;
end

if any(strcmp(varargin, 'no_warnings'))
    doprintwarnings = false;
end

if any(strcmp(varargin, 'cosine_similarity')) % run cosine instead of dot-product
    sim_metric = 'cosine';
end

if any(strcmp(varargin, 'correlation')) || any(strcmp(varargin, 'corr')) % run correlation instead of dot-product
    sim_metric = 'corr';
end

if any(strcmp(varargin, 'binary_overlap')) % run overlap instead of dot-product
    sim_metric = 'overlap';
end

if any(strcmp(varargin, 'dotproduct')) % default. overwrites previous selections
    sim_metric = 'dotproduct';
end

if any(strcmp(varargin, 'posterior_overlap')) % run overlap instead of dot-product
    sim_metric = 'posterior_overlap';
end

if any(strcmp(varargin, 'treat_zero_as_data'))
    treat_zero_as_data = true;
end

if any(strcmp(varargin, 'exclude_zero_mask_values'))
    exclude_zero_mask_values = true;
end

% if docosine && docorr, error('Choose either cosine_similarity or correlation, or no optional inputs for dot product'); end

% ---------------------------------
% Variable types and sizes
% ---------------------------------
pattern_weights = double(pattern_weights); % force double b/c of matlab instabilities
dat = double(dat); % force double b/c of matlab instabilities

[n, k] = size(dat);
[n2, npatt] = size(pattern_weights);

if n ~= n2, error('Number of observations must be equal for data and patterns'); end

similarity_output = NaN .* zeros(k, npatt);

% ---------------------------------
% Missing/excluded values image-wise
% ---------------------------------

if treat_zero_as_data 

    badvals = isnan(dat);

else

    badvals = dat == 0 | isnan(dat);  % Matrix. not used for binary overlap (SG).

end

if exclude_zero_mask_values
    
    badvals_mask = pattern_weights == 0 | isnan(pattern_weights);
    
else

    badvals_mask = isnan(pattern_weights);

end


% ---------------------------------
% Main similarity calculation
% ---------------------------------

% if ~docorr && ~docosine
if strcmp(sim_metric,'dotproduct')
    % dot product. No need to remove missing (zero-valued) voxels because dotproduct is
    % the same either way.
    
    similarity_output = dotproduct(dat, pattern_weights);
    
else
    
    % all other metrics
    for i = 1:npatt
        
        switch sim_metric
            case 'corr'
               
                similarity_output(:, i) = image_correlation(dat, pattern_weights(:, i), badvals, badvals_mask);

            case 'cosine'
                
                similarity_output(:, i) = cosine_similarity(dat, pattern_weights(:, i), badvals);
                
            case 'overlap'
                
                similarity_output(:, i) = overlap_similarity(dat, pattern_weights(:, i));
                
            case 'posterior_overlap'
                
                similarity_output(:, i) = posterior_overlap_similarity(dat, pattern_weights(:, i));
                
            otherwise
                error('Invalid similarity metric.');
        end
        
    end
    
end % pattern sim calculation


% ---------------------------------
% Weight cases based on distance to mean
% ---------------------------------

if any(strcmp(varargin,'weighted')) % if we want to weight data based on similarity to group mean

    mean_dat = mean(dat,2); %sample mean
    distances = squareform(pdist([dat,mean_dat]')); %use other measures?
    weights=1./distances(1:end-1,end);
    weights=weights-min(weights)+1e-12; %min0
    weights=weights./mean(weights); % scale so weights are positive with mean value one

else %
    weights=ones(size(dat,2),1); %just use a weight of one
end


similarity_output = bsxfun(@times,similarity_output,weights);

% ---------------------------------
% Warnings
% ---------------------------------

if doprintwarnings && any(badvals(:)) && ~doignoremissing && ~strcmp(sim_metric,'overlap')
    
    
    for i = 1:npatt
        
        inmask = ~(pattern_weights(:, i) == 0 | isnan(pattern_weights(:, i)));
        
        bad_in_mask = bsxfun(@times, badvals, inmask);
        bad_in_mask = sum(bad_in_mask);  % how many bad values are in mask, across images
        
        if any(bad_in_mask)
            if treat_zero_as_data
                fprintf('Warning: Some data images have zero values in some of the %3.0f valid voxels in weight mask. \nThese are treated as data, affecting similarity values (dot product or cosine sim effectively excludes these voxels but correlation metric is influenced)\n', sum(inmask));
                disp('Number of zero or NaN values within weight mask, by input image:');

            else % zeros are excluded
                fprintf('Warning: Some data images have zero values in some of the %3.0f valid voxels in weight mask. \nThese voxels are excluded from similarity analysis image-wise.\n', sum(inmask));
                disp('Number of zero or NaN values within weight mask, by input image:');
            end

            for j = 1:length(bad_in_mask)
                fprintf('%3.0f ', bad_in_mask(j));
            end
            fprintf('\n');
        end
        
    end  % pattern index
    
end % bad val warnings

end % main function




% ---------------------------------
% *
% Sub-functions
% *
% ---------------------------------


% FUNCTIONS for dot product, cosine_similarity, correlation with missing
% values (voxels) that may vary on a data image-by-image basis.

function pexp = dotproduct(X, weights)

% X = voxels x images data matrix, e.g., dat.dat
% weights is voxels x 1 weight data
% badvals is voxels x images logical matrix of voxels to exclude image-wise

% badvals do not matter
pexp = (weights' * X)';

end % function


function [ab, a, b] = image_norms(dat, pattern_weights, badvals)
% [a, b] = image_norms(dat, pattern_weights, badvals)
% dat = voxels dat images data matrix, e.g., dat
% pattern_weights is voxels dat 1 weight data
% badvals is voxels dat images logical matrix of voxels to exclude image-wise

% Norm for image data.
% All non-zero, non-NaN values valid. zero/nan have no
% effect implictly.

a = nansum(dat .^ 2)' .^ .5;

for i = 1:size(dat, 2)    % Loop because we may have different voxel exclusions in each image
    
    inmask = ~badvals(:, i);
    
    b(i, 1) = nansum(pattern_weights(inmask) .^ 2) .^ .5;  % Norm for pattern_weights, excluding out-of-image pattern_weights image-wise
    
end

ab = a .* b;

end % function


function cossim = cosine_similarity(dat, pattern_weights, badvals)

pexp = dotproduct(dat, pattern_weights);

ab = image_norms(dat, pattern_weights, badvals);

cossim = pexp ./ ab;

% when division by zero (e.g. non-overlapping binary masks) return 0.
% Stephan 2017/3/8
if any(ab==0)
    cossim(isnan(cossim)) = 0;
end

end % function


function r = image_correlation(dat, pattern_weights, badvals, badvals_mask)


for i = 1:size(dat, 2)    % Loop because we may have different voxel exclusions in each image
    
    inmask = ~badvals(:, i) & ~badvals_mask(:, i);
            
    r(i, 1) = corr(pattern_weights(inmask), dat(inmask, i));  % Correlation, excluding out-of-image pattern_weights image-wise
        
   
end


end % function


function r = overlap_similarity(dat, pattern_weights)

% check for binary data
if numel(unique(dat(:)))>2 || sum(unique(dat(:))-[0; 1])~=0 ...
        || numel(unique(pattern_weights(:)))>2 || sum(unique(pattern_weights(:))-[0; 1])~=0
    if numel(unique(round(pattern_weights))) == 2
        pattern_weights = round(pattern_weights); % an easy (but maybe not optimal) solution for the resampled binary mask
    else
        error('Binary overlap similarity needs binary data [0 1] input.');
    end
end

% compute overlap
for i = 1:size(dat, 2)
    
    inmask = isfinite(dat(:,i)) & isfinite(pattern_weights);
    nVox = sum(inmask);
    
    r(i,1) = sum(dat(inmask,i)==1 & pattern_weights(inmask)==1) / nVox; % overlap excluding NaNs
    
end

end % function


function r = posterior_overlap_similarity(dat, pattern_weights)

% check for binary data
if numel(unique(dat(:)))>2 || sum(unique(dat(:))-[0; 1])~=0 ...
        || numel(unique(pattern_weights(:)))>2 || sum(unique(pattern_weights(:))-[0; 1])~=0
    if numel(unique(round(pattern_weights))) == 2
        pattern_weights = round(pattern_weights); % an easy (but maybe not optimal) solution for the resampled binary mask
    else
        error('Binary overlap similarity needs binary data [0 1] input.');
    end
end

% compute posterior overlap
for i = 1:size(dat, 2)
    
    % calculate the space
    inmask = sum(dat,2)~=0;
    % calculate P(pattern | data) using overlaps
    r(i,1) = sum(dat(inmask,i)==1 & pattern_weights(inmask)==1)./sum(dat(inmask,i)==1);
    
end


end

function R = weightedcorrs(Y, w)
%
%   WEIGHTEDCORRS returns a symmetric matrix R of weighted correlation
%   coefficients calculated from an input T-by-N matrix Y whose rows are
%   observations and whose columns are variables and an input T-by-1 vector
%   w of weights for the observations. This function may be a valid
%   alternative to CORRCOEF if observations are not all equally relevant
%   and need to be weighted according to some theoretical hypothesis or
%   knowledge.
%
%   R = WEIGHTEDCORRS(Y, w) returns a positive semidefinite matrix R,
%   i.e. all its eigenvalues are non-negative (see Example 1).

%
% % ======================================================================
%
%   See also CORRCOEF, COV, STD, MEAN.
%
% % ======================================================================
%
%-*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-*%
%                                                                                               %
%            Author: Liber Eleutherios                                             %
%            E-Mail: libereleutherios@gmail.com                             %
%            Date: 23 July 2008                                                      %
%            Updated: 6 June 2012                                                 %
%                                                                                               %
%-*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-*%
%
% % ======================================================================
%
% Check input
ctrl = isvector(w) & isreal(w) & ~any(isnan(w)) & ~any(isinf(w)) & all(w > 0);
if ctrl
    w = w(:) / sum(w);                                                          % w is column vector
else
    error('Check w: it needs be a vector of real positive numbers with no infinite or nan values!')
end
ctrl = isreal(Y) & ~any(isnan(Y)) & ~any(isinf(Y)) & (size(size(Y), 2) == 2);
if ~ctrl
    error('Check Y: it needs be a 2D matrix of real numbers with no infinite or nan values!')
end
ctrl = length(w) == size(Y, 1);
if ~ctrl
    error('size(Y, 1) has to be equal to length(w)!')
end
[T, N] = size(Y);                                                             % T: number of observations; N: number of variables
temp = Y - repmat(w' * Y, T, 1);                                              % Remove mean (which is, also, weighted)
temp = temp' * (temp .* repmat(w, 1, N));                                     % Covariance Matrix (which is weighted)
temp = 0.5 * (temp + temp');                                                  % Must be exactly symmetric
R = diag(temp);                                                               % Variances
R = temp ./ sqrt(R * R');

end % function


