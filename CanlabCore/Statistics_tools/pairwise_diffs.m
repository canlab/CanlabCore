function d = pairwise_diffs(x, varargin)
% Pairwise operations on columns of a matrix x, for each row
%
% :Usage:
% ::
%
%     d = pairwise_diffs(x, [function handle])
%
% :Inputs:
%
%   **x:**
%        an n x k matrix
%
%   Optional: a function handle for the operation to perform on pairwise
%   elements of x(i, :)
%
% :Outputs:
%
%   if x is n x k, d is n x (k*(k-1)/2)
%
%   columns of x are arranged this way, e.g., with k = 5:
%
%   [x(1, 1) - x(1, 4)  x(1, 1) - x(1, 5)  x(1, 2) - x(1, 3) ...]
%
%   squareform(d(1, :)) is a matrix of pairwise diffs for d
%
% :Examples:
% ::
%
%    x = magic(5)                  % generate data
%    d = pairwise_diffs(x);        % d = pairwise differences
%    row1 = squareform(d(1, :));   % square matrix of pairwise diffs for row 1
%
%    d = pairwise_diffs(x, @(a, b) a + b);  % return pairwise sum instead
%
% ..
%    Tor Wager, Nov 2012
% ..

fhan = @(a, b) a - b; % default: difference operator

if length(varargin) > 0, fhan = varargin{1}; end

for i = 1:size(x, 1)

    d(i, :) = get_result(x(i, :), fhan);
    
end


end % function


% Subfunctions

function d = get_result(x, fhan)

    n = length(x);
    indx = 1;
    
    for i = 1:n
    
        for j = i+1:n
           
            d(1, indx) = fhan(x(i), x(j));
            indx = indx + 1;
            
        end
        
    end
    

end

