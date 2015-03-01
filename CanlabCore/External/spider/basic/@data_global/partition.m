function [d1, d2, sel] = partition(d, sel, p)

% PARTITION - partition a DATA_GLOBAL object randomly or deterministically
% 
% D1 = PARTITION(D, SEL)
% [D1, D2] = PARTITION(D, SEL)
% [D1, D2, SEL] = PARTITION(D, 'rand', P)
% 
% The same as the DATA/PARTITION method, except that it operates on
% DATA_GLOBAL objects, manipulating D.index, and D.myX and D.myY if
% present.

n = length(d.index);

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

d1 = d;
d1.index(~sel) = [];
if ~isempty(d1.myX), d1.myX(~sel, :) = []; end
if ~isempty(d1.myY), d1.myY(~sel, :) = []; end

if nargout < 2, return, end

d2 = d;
d2.index(sel) = [];
if ~isempty(d2.myX), d2.myX(sel, :) = []; end
if ~isempty(d2.myY), d2.myY(sel, :) = []; end
