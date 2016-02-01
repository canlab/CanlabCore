function [yhat,yhat_sorted,err,goodness,mymad] = monotonic_regression(x,y,doplot)
% See lsqisotonic.m
%
% ..
%    tor wager, nov. 06
%
%    mymad is median absolute deviation from isotonic regression line
%    goodness was mad(yhat) / mymad : ratio of mean values to error
%    is now % variance explained by fit (r-square)
% ..

if nargin < 3, doplot = 0; end

if doplot
    figure('Color','w');plot(x,y,'ko','MarkerFaceColor','k')
end

n = numel(x);

[xyord,ord] = sortrows([x(:) y(:)]); 
iord(ord) = 1:n;
xyord = double(xyord);
% Initialize fitted values to the given values.
yhat = xyord(:,2);
w = ones(size(yhat));

if(n <= intmax('uint16'))
    dtconv = @uint16;
elseif(n <= intmax('uint32'))
    dtconv = @uint32;
else
    dtconv = @double;
end

block = dtconv(1:n);


while true
    % If all blocks are monotonic, then we're done.
    diffs = diff(yhat);
    if all(diffs >= 0), break; end

    % Otherwise, merge blocks of non-increasing fitted values, and set the
    % fitted value within each block equal to a constant, the weighted mean
    % of values in that block.
    idx = dtconv(cumsum([1; (diffs>0)]));
    sumyhat = accumarray(idx,w.*yhat);
    w = accumarray(idx,w);
    yhat = sumyhat ./ w;
    block = idx(block);
end


yhat_sorted = yhat(block);

if doplot
    hold on; plot(xyord(:,1),yhat_sorted,'ro-','LineWidth',2);
end

yhat = reshape(yhat_sorted(iord), size(y));
err = y - yhat;
if nargout > 4, mymad = mad(err); end
%goodness = mad(yhat) ./ mymad;
vv = var(y);
goodness = 1 - (var(err) ./ vv);

return
