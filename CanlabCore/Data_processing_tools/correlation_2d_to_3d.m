function C = correlation_2d_to_3d(Cmatrix)
% Turn a 2-d flat matrix of correlations (in squareform.m format) into an k x k x n matrix of correlations 
% 
% C = correlation_2d_to_3d(Cmatrix)
%
% C = [k x k x n],  each [k x k] matrix is inter-variable correlation
% matrix, one correlation matrix for each of n replicates.

% get data for prediction
%-----------------------------------

[n, k] = size(Cmatrix);

E = eye(k);

Cmatrix = zeros(k, k, n);

for i = 1:n
    
    m = squareform(Cmatrix(i, :));
    m = m + m' + E;
    
    C(:, :, i) = m;
end

end
