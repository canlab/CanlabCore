function hoxlist = hoxsetup(numhox,hoxrange,sizeGenerations)
%
% Tor Wager, 12/29/01

if numhox ~= size(hoxrange,1), error('num hox genes does not equal number of rows of hoxrange.'),end

hoxlist = rand(numhox, sizeGenerations);

for i = 1:size(hoxlist,1)
    hoxlist(i,:) = hoxlist(i,:) * (hoxrange(i,2) - hoxrange(i,1)) + hoxrange(i,1);
end

return