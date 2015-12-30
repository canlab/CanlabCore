function [x,V,xc] = smooth_timeseries(x,perc)
% Make exponential smoothing function with perc proportion of data points
% (0 < perc < 1)
% OR specify length directly, as % of points to 0 weight
%
% :Usage:
% ::
%
%     [x,V] = smooth_timeseries(x,perc)
%
% apply smoothing filter V to data (x), V * x
% (works just as well for a matrix of column vectors)
% You could also apply it to a model matrix X

if size(x,1) == 1, x = x';,end  % transpose if necessary

if perc > 1, len = perc;,
else
    len = ceil(size(x,1) .* perc);
end
xc = 1:len;
xc = 1 ./ (1 + xc.^.5);
xc = xc - min(xc);
%xc = xc ./ sum(xc); irrelevant

V = toeplitz(pad(xc',size(x,1)-length(xc)));
V = V ./ repmat(sum(V),size(V,1),1);
V = V';

x = V * x;

return
