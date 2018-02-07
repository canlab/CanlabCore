%% lighten underlay volume - for white background
function Z = lighten_underlay_edges(Z, maxpercent)

Zi = Z(:);

mindat = min(Zi);
maxdat = max(Zi);

Z(Z == mindat) = maxdat;

% higher ending values of k will be more 'softening'; lower =
% more dark edges

myprctiles = prctile(Zi, [0:1:maxpercent]); % speed up!!

for k = 1:maxpercent % for each of the darkest/lowest-value percentiles
    
    w = 1 - k/100; % weights, for wtd average of orig and inverted, 1 is all inverted
    % (bilinear weighting)
    
    %isdark = Zi < prctile(Zi, k) & Zi >= prctile(Zi, k - 1);
    isdark = Zi < myprctiles(k + 1) & Zi >= myprctiles(k);
    
    Z(isdark) = maxdat;
    
    %Z(isdark) = w .* (maxdat) + (1 - w) .* Z(isdark);
    
end % the weight value loop

end