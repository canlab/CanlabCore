function s = stress(d,dhat,varargin)
% s = stress(d,dhat,[k])
%
% s     = stress
% d     = observed distances
% dhat  = implied distances
% [k]   = optional # dimensions to save
%
% if 3rd argument is entered, dhat is assumed to be a matrix of stimulus
% coordinates (e.g., output of cmdscale), and the 3rd argument is the
% number of columns to use for computing implied distances.
%
% From Shepard, 1980.
% by Tor Wager, 8/1/04

if length(varargin) > 0
    dhat = squareform(pdist(dhat(:,1:varargin{1})));
end


% numerator
dif = (d - dhat) .^ 2;
dif = sum(dif(:));

% denominator
den = sum(sum(d .* d));

s = dif / den;

return