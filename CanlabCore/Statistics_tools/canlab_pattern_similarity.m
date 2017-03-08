function similarity_output = canlab_pattern_similarity(dat, pattern_weights, varargin)
% Calculate similarity between each column in a data matrix dat and a vector of pattern weights
%
% - Similarity options: dot product, cosine similarity, and correlation
% - Columns are often images, e.g., from fmri_data.dat for fmri_data objects
% - weights are pattern weights from one or more 'signature' patterns (each pattern is a column)
%
% - Assumes dat and pattern_weights matrices have equal rows (voxels) and include valid voxels
% - Removes empty column data column-wise (image-wise) in case some images have uneven voxel coverage
%
% - Used in apply_mask, extract_roi_averages, image_similarity_plot***
%
% similarity_output = canlab_pattern_similarity(dat, pattern_weights)
%
%
% :Optional Inputs:
%
%   **cosine_similarity**
%       Use cosine similarity metric for pattern expression instead of dot product
%
%   **correlation**
%       Use correlation metric for pattern expression instead of dot product
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
% Pattern pattern_weights can have zeros, which may be valid values in voxels,
% i.e., with binary masks
% Images with values of zero or NaN are considered out-of-mask, as they
% are not valid values. These should be excluded from both image and
% mask when calculating similarity.
% Thus, there is an asymmetry between pattern mask and image data
% in considering which voxels to use.
% Otherwise, all dot product and similarity metrics are standard.

% We also need to calculate bad values on an image-by-image basis, not
% relying on remove_empty to exclude voxels with ineligible values
% across the entire set of images.
%
%  Created by Tor Wager - 3/7/17


% ---------------------------------
% Defaults and optional inputs
% ---------------------------------

docorr = false;     % run correlation instead of dot-product for pattern expression
docosine = false;   % run cosine sim instead
doignoremissing = false; % ignore warnings for missing voxels

if any(strcmp(varargin, 'ignore_missing'))
    doignoremissing = true;
end

if any(strcmp(varargin, 'cosine_similarity'))
    docosine = true;
end

if any(strcmp(varargin, 'correlation')) % run correlation instead of dot-product
    docorr = true;
end

if docosine && docorr, error('Choose either cosine_similarity or correlation, or no optional inputs for dot product'); end

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

badvals = dat == 0 | isnan(dat);  % Matrix

% ---------------------------------
% Main similarity calculation
% ---------------------------------

for i = 1:npatt
    
    if docorr
        
        similarity_output(:, i) = image_correlation(dat, pattern_weights(:, i), badvals);
        
    elseif docosine
        
        similarity_output(:, i) = cosine_similarity(dat, pattern_weights(:, i), badvals);
        
    else
        
        similarity_output(:, i) = dotproduct(dat, pattern_weights(:, i), badvals);
        
    end
    
    % ---------------------------------
    % Warnings
    % ---------------------------------
    
    if any(badvals(:)) && ~doignoremissing
        
        disp('WARNING!!! SOME SUBJECTS HAVE ZERO VALUES WITHIN IMAGE MASK.');
        disp('This could artifactually influence their scores if these 0 values are in the weight mask. ');
        disp('They will be excluded from similarity analysis image-wise.');
        
        wh = find(any(badvals)); % which images
        
        inmask = ~(pattern_weights(:, i) == 0 | isnan(pattern_weights(:, i)));
        
        fprintf('Total voxels in weight mask: %3.0f\n', sum(inmask));
        disp('Test images with bad values:');
        for i = 1:length(wh)
            fprintf('Test image %3.0f: %3.0f zero values\n', wh(i), sum(badvals(:, wh(i))));
        end
        
    end
    
    
end  % pattern index


end % main function




% ---------------------------------
% *
% Sub-functions
% *
% ---------------------------------


% FUNCTIONS for dot product, cosine_similarity, correlation with missing
% values (voxels) that may vary on a data image-by-image basis.

function pexp = dotproduct(X, weights, badvals)

% X = voxels x images data matrix, e.g., dat.dat
% weights is voxels x 1 weight data
% badvals is voxels x images logical matrix of voxels to exclude image-wise

for i = 1:size(X, 2)    % Loop because we may have different voxel exclusions in each image
    
    inmask = ~badvals(:, i);
    
    pexp(i, 1) = X(inmask, i)' * weights(inmask, 1);
    
end

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

end


function cossim = cosine_similarity(dat, pattern_weights, badvals)

pexp = dotproduct(dat, pattern_weights, badvals);

ab = image_norms(dat, pattern_weights, badvals);

cossim = pexp ./ ab;


end % function


function r = image_correlation(dat, pattern_weights, badvals)

for i = 1:size(dat, 2)    % Loop because we may have different voxel exclusions in each image
    
    inmask = ~badvals(:, i);
    
    r(i, 1) = corr(pattern_weights(inmask), dat(inmask, i));  % Correlation, excluding out-of-image pattern_weights image-wise
    
end


end


