function M=divide_sparse(M, by)
%DIVIDE_SPARSE Divide each row in the sparse matrix M by the number in the
% corresponding row of "by". In practice, this is much faster for large
% arrays than using normal MATLAB syntax.
%
% M = DIVIDE_SPARSE(M, by)
%
% Parameters
% ----------
% M: sparse matrix of size (m1, m2)
% 
% by: array of size (m1, 1)
% 
% Returns
% -------
% M: sparse matrix of size (m1, m2)
%     The result of diving the (i, j)-th entry of M by the i-th component of
%     "by".
%
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

[rows, cols]=find(M);
R=size(rows);

for r=1:R
    M(rows(r), cols(r))=M(rows(r), cols(r))/by(rows(r));
end

end