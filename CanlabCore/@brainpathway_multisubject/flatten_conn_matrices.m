% Given a 3d matrix, n_nodes x n_nodes * subjects, will return a 2D matrix, subjects x n_node-pairs
% by default, will use bs.connectivity.regions.r
% if a matrix is passed in, will use that instead
%
% To revert a 1 x node-pairs matrix back to square form, see square_form.m in CanlabCore
%
% Yoni Ashar, March 2020
%
function flat_conn_mat = flatten_conn_matrices(bs, varargin)

    mat3d = bs.connectivity.regions.r;
    if nargin > 1 & ismatrix(varargin{1})
        mat3d = varargin{1};
    end

    n_subjects = size(mat3d, 3);
    n_nodes = size(mat3d, 2);

    % there is probably a faster way to do this
    for i=1:n_subjects
        tmp = tril(mat3d(:,:,i), -1);
        flat_conn_mat(i,:) = tmp(tril(true(n_nodes), -1)); % could also use squareform(), but this requires 0's along the diagonal, so would need to set that first. no clear benefit to using that function so will leave as is. Yoni March 2020
    end
end