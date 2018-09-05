function [D2, D2_expected, pval, wh_outlier_uncorr, wh_outlier_corr] = mahal(obj, varargin)
% Mahalanobis distance for each image in a set compared to others in the set
% Lower p-values reflect higher outlier status for an image.
%
% [D2, D2_expected, pval, wh_outlier_uncorr, wh_outlier_corr] = mahal(obj, varargin)
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

