function [D2, D2_expected, pval, wh_outlier_uncorr, wh_outlier_corr] = mahal(obj, varargin)
% Mahalanobis distance for each image in a set compared to others in the set
% Lower p-values reflect higher outlier status for an image.
%
% [D2, D2_expected, pval, wh_outlier_uncorr, wh_outlier_corr] = mahal(obj, varargin)
%
% Mahalanobis distance, which is a measure of how different each image is from the rest of the images. 
% It's based on the squared distance across images, just like the typical least-squares solution we 
% use to fit regression models (the basis for all General Linear Models). But the Mahalanobis distance
% provides distances along the principal axes of variation based on the covariance across images. 
% Here is an intuition for why this is important. A difference of a given magnitude doesn't mean the 
% same thing for all voxels, or for all variables in a dataset in general.  If I'm classifying people 
% in terms of their height in inches weight in lbs, and the number of fingers they have, +1 inch in 
% height isn't so surprising, but +1 finger is. And if I transform height to feet, +1 foot is also 
% very surprising. Rather than differences in raw units, it often makes sense to measure deviations 
% in units of standard devations, considering the intrinsic variablilty in each variable. 
% Furthermore, variables can move together (i.e., covary). If someone is +1 sd in height, they're 
% also likely to be heavier than average as well, so I wouldn't want to sum the deviation in height 
% and deivation in weight without considering their covariance and discounting it. If height and 
% weight covary positively, +1 height and -1 weight would be much more surprising, because this pattern 
% goes against the natural covariance in these measures. Mahalanobis distance takes care of all of that, 
% by calculating deviations along the principal axes of covariation. 
%
% A note on outlier identification:
% wh_outlier_uncorr identifies images that are outside the 95% confidence region of the cloud of 
% images (think of each image as a point) in multidimensional space. This is an indicator vector 
% (values of 1 and 0) for which images are outliers. The plot shows us the distances for each case (image)
% and the relationship between the observed and expected distance. A positive deviation from the line means 
% the image is more unusual than expected.
% 
% Optional inputs:
% 'corr' : use correlation matrix of images instead of covariance. 
% Insensitive to differences in scale and mean of data, and thus more
% sensitive to changes in pattern across image.
% 
% Outputs:
%
% wh_outlier_uncorr
%   Logical vectors of which images have mahalanobis distances with p < .05
%   and values greater than the median
%
% wh_outlier_corr
%   Logical vectors of which images have mahalanobis distances with p < .05
%   corrected (Bonferroni), and values greater than the median
%
% Examples:
% ----------------------------------------------------------------------
% [ds, expectedds, p, wh_outlier_uncorr, wh_outlier_corr] = mahal(fmridat, 'noplot');
% 
% Y = ds - expectedds;
% wh = p < (.05 ./ length(p));  % Outliers after Bonferroni correction
% 
% plot(Y);
% plot(find(wh), Y(wh), 'ro', 'MarkerSize', 6);
%
% Tor Wager
%
% 

% Programmers' notes: tor - 4/6/2018, changed num of components retained to
% avoid empty components, evaluated as better performance against cov
% plots. Also added 'corr' mode to work on correlation rather than cov
% matrix [optional]

doplot = true;
if any(strcmp(varargin, 'noplot')), doplot = false; end

if any(strcmp(varargin, 'corr'))
    
    mydata = double(obj.dat);
    mydata = zscore(mydata)';
    
else
    
    mydata = double(obj.dat');
    
end


[coeff, score, latent, tsquared, explained] = pca(mydata, 'Economy', true);

% Taking up to 90% of var explained was adding too much noise - mahal not
% good with empty dimensions.  Reduce further: Each pc retained must
% explain > 1 observations's worth if all obs were independent, ie., 100/n
% below. 
%wh = find(cumsum(explained) <= 90);

n = length(explained);
wh = explained > 100/n;


d = gradient(latent);
rapidchange = d ./ d(1) > .05; % gradint is at least 5% of initial drop 

wh = wh & rapidchange;

fprintf('Retained %d components for mahalanobis distance\n', sum(wh));

% matlab's version: same, but no p-values, etc.
% D2 = mahal(score(:, wh), score(:, wh));

if doplot
    
    [ds, S, pval] = multivar_dist(score(:, wh));
    
else
    
    [ds, S, pval] = multivar_dist(score(:, wh), 'noplot');
    
end

D2 = ds(:, 1);
D2_expected = ds(:, 3);

%figure; plot(ds(:, 3), ds(:, 1), 'ro');

wh_outlier_uncorr = pval < .05 & D2 > median(D2);

wh_outlier_corr = pval < (.05 ./ length(pval)) & D2 > median(D2);  % Outliers after Bonferroni correction


fprintf('\n')
fprintf('Potential outliers based on mahalanobis distance:\n')
fprintf('Bonferroni corrected: %d images\t\tCases ', sum(wh_outlier_corr));
fprintf('%d ', find(wh_outlier_corr))
fprintf('\n')
fprintf('Uncorrected: %d images\t\tCases ', sum(wh_outlier_uncorr));
fprintf('%d ', find(wh_outlier_uncorr))
fprintf('\n')
fprintf('\n')

end

