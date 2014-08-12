function [my_ste,t,n,p,m] = ste(dat)
% function [my_ste,t,n_in_column,p,mean] = ste(dat)
%
% tor wager 
%
% standard error of the mean
% and t-values, columnwise
% and p-values, if asked for
%
% omits NaN values row-wise within each column
%  (edit Nov 30, 2006 by tor)
%
% does NOT use n - 1
% matches matlab t-test function

n = sum(~isnan(dat),1);
%if ~isempty(whom), dat(whom,:) = [];, end

m = nanmean(dat);
my_ste = nanstd(dat) ./ (n.^.5);

if nargout > 1
    t = m ./ my_ste;
end

if nargout > 3
    p = 2 .* (1 - tcdf(abs(t), n - 1));
end

return