function Cmatrix = correlation_3d_to_2d(C)
% Turn an k x k x n matrix of correlations into a 2-d flat matrix
% 
% Cmatrix = correlation_3d_to_2d(C)
%
% C = [k x k x n],  each [k x k] matrix is inter-variable correlation
% matrix, one correlation matrix for each of n replicates.

% get data for prediction
%-----------------------------------

[j, k, n] = size(C);

if j ~= k, error('Not a series of correlation matrices'); end

E = eye(j);

nout = j * (j - 1) ./ 2;
Cmatrix = zeros(n, nout);

for i = 1:n
    
    Cmatrix(i, :) = squareform(C(:, :, i) - E);
    
end

end
