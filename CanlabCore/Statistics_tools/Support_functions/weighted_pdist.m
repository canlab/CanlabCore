function d = weighted_pdist(x,w)
% d = weighted_pdist(x,w)
%
% x is matrix of locations, rows are objects, cols are dimensions, vals are
% coordinates.
%
% w is a row vector of weights for dims (cols) of x
%   weights are applied to squared distances
%   sqrt(  sum over dims (w * (xi - xj).^2 )  )
%
% d is a square matrix of Euclidean distances btwn points
% tor wager

%w = w .^ 2;  this gives weight pattern as if weights are applied as XW

for i = 1:size(x,1)
    for j = i+1:size(x,1)
        
        % squared distance
        tmp = (x(i,:) - x(j,:)) .^ 2;   % squared distances
        
        tmp = tmp * w';                 % weighted sum
        
        d(i,j) = tmp .^ .5;
        
    end
end

d(end+1,:) = 0;
d = d + d';

