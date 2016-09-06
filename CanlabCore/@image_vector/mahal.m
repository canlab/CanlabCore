function [D2, D2_expected, pval] = mahal(obj, varargin)

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

