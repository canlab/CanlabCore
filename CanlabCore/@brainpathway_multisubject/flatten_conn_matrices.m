% Given a 3d matrix, n_nodes x n_nodes * subjects, will return a 2D matrix, subjects x n_node-pairs
%
% flat_conn_mat = flatten_conn_matrices(bs, ['replacenans'])
%
% by default, will use bs.connectivity.regions.r
% if a matrix is passed in, will use that instead
%
% To revert a 1 x node-pairs matrix back to square form, see square_form.m in CanlabCore
%

% Programmers' Notes:
% Yoni Ashar, March 2020
% Tor Wager, Dec 2020 -> edited to add support for NaN imputation, and
%                           speed up by pre-defining output matrix

function flat_conn_mat = flatten_conn_matrices(bs, varargin)

    doreplacenans = false;
    if any(strcmp(varargin, 'replacenans'))
        doreplacenans = true;
    end
    
    mat3d = bs.connectivity.regions.r;
    
    if nargin > 1 && ismatrix(varargin{1}) && ~ischar(varargin{1})
        
        mat3d = varargin{1};
        
    end

    n_subjects = size(mat3d, 3);
    n_nodes = size(mat3d, 2);

    % pre-load to make faster
    flat_conn_mat = zeros(n_subjects, n_nodes * (n_nodes - 1) ./ 2);
    
    for i=1:n_subjects
        
        tmp = tril(mat3d(:,:,i), -1);
        flat_conn_mat(i,:) = tmp(tril(true(n_nodes), -1)); % could also use squareform(), but this requires 0's along the diagonal, so would need to set that first. no clear benefit to using that function so will leave as is. Yoni March 2020
    
    end
    
    if doreplacenans
        % Replace NaNs with subject-wise means
        
        % Find NaNs and replace with subject mean (so regions will have min=zero influence on isc)
        fprintf('Imputing missing (NaN) values. ')
        
        whnan = isnan(flat_conn_mat);
        subjmean = nanmean(flat_conn_mat');
        
        for i = 1:size(flat_conn_mat, 1)
            flat_conn_mat(i, whnan(i, :)) = subjmean(i);
        end
        
    end
    
end % main function