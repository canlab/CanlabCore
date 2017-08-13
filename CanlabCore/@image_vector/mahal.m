function [D2, D2_expected, pval] = mahal(obj, varargin)
% [D2, D2_expected, pval] = mahal(obj, varargin)
%
% Mahalanobis distance for each image in a set compared to others in the set
% Lower p-values reflect higher outlier status for an image.
%
% Examples:
% ----------------------------------------------------------------------
% [ds, expectedds, p] = mahal(fmridat, 'noplot');
% 
% Y = ds - expectedds;
% wh = p < (.05 ./ length(p));  % Outliers after Bonferroni correction
% 
% plot(Y);
% plot(find(wh), Y(wh), 'ro', 'MarkerSize', 6);
%
% Tor Wager

doplot = true;
if any(strcmp(varargin, 'noplot')), doplot = false; end

[coeff, score, latent, tsquared, explained] = pca(obj.dat', 'Economy', true);
wh = find(cumsum(explained) <= 90);

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

end

