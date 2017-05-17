function [d1, d2, sel] = partition(d, sel, p)

% PARTITION - partition a DATA object randomly or deterministically
% 
% D1 = PARTITION(D, SEL)
% [D1, D2] = PARTITION(D, SEL)
% 
% D is a DATA object. SEL is either a list of indices to data points in D,
% or a LOGICAL vector with one element per data point.
% 
% D1 returns a DATA object containing the subset of data points from D
% selected by SEL. D2, if requested, contains the remaining data points. 
% 
% [D1, D2, SEL] = PARTITION(D, 'rand', P)
% 
% Using the 'rand' option, the scalar argument P is either a probability
% (points are assigned to D1 with independent probability P) or an integer
% number of points (a random selection of exactly P points is drawn without
% replacement for assignment to D1). The default value of P is 0.5.
% 
% The third output argument, SEL, is always an Nx1 logical vector
% containing 1 for points that were selected for D1, and 0 for D2.

X = get_x(d);
Y = get_y(d);
n = size(X, 1);
sel = sel(:);
if ischar(sel)
    if ~strcmp(lower(sel'), 'rand'), error(['unknown option ''' sel' '''']), end
    if nargin < 3, p = 0.5; end
    if ~isnumeric(p) | prod(size(p)) ~= 1, error('P must be a numeric scalar'), end
    p = full(real(double(p)));
    if p < 0
        error('P may not be negative')
    elseif p > n
        error('P may not be larger than the number of data points')
    elseif p < 1
        sel = (rand(n, 1) < p);
    elseif p ~= round(p)
        error('P must either be a probability, or a positive integer')
    else
        sel = randperm(n)';
        sel = sel(1:p);
    end
end
if islogical(sel)
    if length(sel) ~= n
        error('logical array must have one element per data point')
    end    
elseif ~isreal(sel) | ~isnumeric(sel)
    error('SEL should be a logical array or an array of valid indices')
else
    sel = full(double(sel));
    if any(sel < 1 | sel > n | sel ~= round(sel))
        error('invalid/out of range indices')
    end
    z = logical(zeros(n,1));
    z(sel) = 1;
    sel = z;
end
if isempty(Y)
    d1 = data(X(sel, :));
else
    d1 = data(X(sel, :), Y(sel, :));
end
if nargout < 2, return, end
if isempty(Y)
    d2 = data(X(~sel, :));
else
    d2 = data(X(~sel, :), Y(~sel, :));
end
