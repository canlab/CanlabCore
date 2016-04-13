%  nonmonotonic_regression_plot(x,y)
% see lsqisotonic.m
function nonmonotonic_regression_plot(x,y)

figure('Color','w');plot(x,y,'ko','MarkerFaceColor','k')

n = numel(x);

[xyord,ord] = sortrows([x(:) y(:)]); iord(ord) = 1:n;
xyord = double(xyord);
% Initialize fitted values to the given values.
yhat = xyord(:,2);
w = ones(size(yhat));

block = 1:n;


while true
    % If all blocks are monotonic, then we're done.
    diffs = diff(yhat);
    if all(diffs >= 0), break; end

    % Otherwise, merge blocks of non-increasing fitted values, and set the
    % fitted value within each block equal to a constant, the weighted mean
    % of values in that block.
    idx = cumsum([1; (diffs>0)]);
    sumyhat = accumarray(idx,w.*yhat);
    w = accumarray(idx,w);
    yhat = sumyhat ./ w;
    block = idx(block);
end


yhat = yhat(block);

hold on; plot(xyord(:,1),yhat,'ro-','LineWidth',2);
yhat = reshape(yhat(iord), size(y));
