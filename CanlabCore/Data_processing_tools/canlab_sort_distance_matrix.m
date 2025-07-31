function [D_sorted, perm] = canlab_sort_distance_matrix(D, varargin)
% canlab_sort_distance_matrix: Reorder distance or correlation matrix so similar items are adjacent
%
% :Usage:
% ::
%     [D_sorted, perm] = canlab_sort_distance_matrix(D, 'method', 'average', 'correlation_matrix', true)
%
% :Inputs:
%   **D:**
%       A [n x n] symmetric distance matrix (or correlation matrix if 'correlation_matrix' = true)
%
% :Optional Inputs:
%   **'method'** (string):
%       Linkage method for hierarchical clustering. Default = 'average'. Options include 'single', 'complete', etc.
%
%   **'correlation_matrix'** (logical):
%       If true, input is a correlation matrix (R). It will be converted to distance via:
%       D = 1 - ((R + 1) / 2). The final output will be re-converted to a correlation matrix.
%       Default = false.
%
%   **'verbose'** (logical):
%       Print info to command window. Default = true.
%
% :Outputs:
%   **D_sorted:**
%       Reordered distance or correlation matrix
%
%   **perm:**
%       Permutation vector: D_sorted = D(perm, perm)
%
% :Examples:
% ::
%     R = corr(randn(100, 20));
%     [R_sorted, perm] = canlab_sort_distance_matrix(R, 'correlation_matrix', true);
%
% :See also:
%   squareform, linkage, optimalleaforder, corr
%


% --------------------------
% Parse inputs
% --------------------------
p = inputParser;
p.addRequired('D', @(x) validateattributes(x, {'numeric'}, {'square', '2d'}));
p.addParameter('method', 'average', @(x) ischar(x) || isstring(x));
p.addParameter('correlation_matrix', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('verbose', true, @(x) islogical(x) || isnumeric(x));
p.parse(D, varargin{:});

method = lower(string(p.Results.method));
is_corr = logical(p.Results.correlation_matrix);
verbose = logical(p.Results.verbose);

% --------------------------
% Convert correlation to distance if needed
% --------------------------
if is_corr
    if verbose, fprintf('Converting correlation matrix to distance matrix...\n'); end
    D = 1 - ((D + 1) ./ 2);  % Map R [-1,1] to D [0,1]
end

% --------------------------
% Convert to condensed vector form
% --------------------------
dVec = squareform(D);

% --------------------------
% Hierarchical clustering
% --------------------------
Z = linkage(dVec, method);
leafOrder = optimalleaforder(Z, dVec);

% --------------------------
% Reorder the matrix
% --------------------------
D_sorted = D(leafOrder, leafOrder);
perm = leafOrder;

% --------------------------
% Convert back to correlation matrix if needed
% --------------------------
if is_corr
    if verbose, fprintf('Converting sorted distance back to correlation...\n'); end
    D_sorted = 2 * (1 - D_sorted) - 1;
end

end