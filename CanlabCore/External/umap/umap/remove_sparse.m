function M=remove_sparse(M, op)
%REMOVE_SPARSE Set to zero the entries of M for which op(M) = 1.
%
% M = REMOVE_SPARSE(M, op)
%
% Parameters
% ----------
% M: sparse matrix of size (m1, m2)
% 
% op: a function accepting doubles as input and outputting a boolean.
% 
% Returns
% -------
% M: sparse matrix of size (m1, m2)
%     The result of setting to zero the entries of M for which op(M) = 1.
%
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

idxs=find(M);
logicalIdxs=feval(op, M(idxs));
removeIdxs=idxs(logicalIdxs);
M(removeIdxs)=0;
end