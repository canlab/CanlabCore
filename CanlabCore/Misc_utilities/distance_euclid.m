function d = distance_euclid(u, v)
% Euclidean distance between u and v
%
% :Usage:
% ::
%
%     d = distance_euclid(u, v)
%
% :Inputs:
%
%   **u** and **v:**
%        should be column vectors or matrices in k-dim space, where k is
%        the number of columns of u and v.
%
%        if u is a single row, replicates u to size(v,1)

if size(u,1) < size(v,1)
    u = repmat(u,size(v,1),1);
    if size(u,1) < size(v,1), error('Inputs must have same number of rows, or one row for 1st input!'),end
end

delt = u - v;
d = (sum(delt .^ 2,2)) .^.5;

return

