% Given a vector of size 1 x n node-pairs, revert back it back to a square
% form (n nodes x n nodes).
%
% Input:
%   vector
%   size of square matrix (just one number)
%   optional value to set for diagonal
%
% Works well to complement brainpathways_multisubject.flatten_conn_matrices
%
% Yoni Ashar, March 2020
%
function square = square_form(vec, n_nodes, varargin)

square = tril(ones(n_nodes), -1); % sets indices for assignment in next line
square(square>0) = vec; % assign in values in correct (columnwise) order
square = square + square'; % flip to be mirrored on upper tri too

% set diagonal values
if ~isempty(varargin) && isnumeric(varargin{1})
    square(logical(eye(size(square)))) = varargin{1}; % set the diagonal p values to be 1
end

end