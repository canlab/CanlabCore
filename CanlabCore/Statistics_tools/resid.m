function [r,X] = resid(X,y, varargin)
% [r,X] = resid(X,y, [add mean back in flag])
%
% tor wager
% residuals from model fit
%
% adds intercept, if missing; adds mean response for each variable if
% requested
% Last edited: 6/2013

add_int = 0;
if ~isempty(varargin) > 0 && varargin{1}
    add_int = 1;
end

% find intercept
wint = all(X == repmat(X(1,:), size(X,1), 1));
if ~any(wint)
    X(:,end+1) = 1; 
    wint = zeros(1, size(X, 2));
    wint(end) = 1;
end

y = double(y);

% break up for efficiency - otherwise Matlab chokes with large
    % matrices...(why?)
    px = pinv(X);
    pxy = px * y;
    xpxy = X * pxy;

    
r = y - xpxy; % X * pinv(X) * y;

if add_int
    m = mean(y);
    ym = repmat(m, size(X, 1), 1);
    r = r + ym;
end

end
