function stats = signtest_matrix(dat)
% :Usage:
% ::
%
%     stats = signtest_matrix(dat)
%
% This is a matric-ized version of Matlab 2010's signtest.m
% it returns identical P-values to Matlab's function
%
% ..
%    copyright 2011, tor wager
% ..

n = size(dat, 1); % n: how many obs in sign test
n_omit = sum(dat == 0 | isnan(dat));
n = (n - n_omit)'; 

isvalid = n > 0;

dat(:, ~isvalid) = [];
n(~isvalid) = [];

npos = sum(dat > 0);
nneg = sum(dat < 0);
sgn = min(npos, nneg)';
mydirection = double(npos > nneg)';
mydirection(mydirection == 0) = -1;

if max(n) < 100
    method = 'exact';
else
    method = 'approximate';
end


if isequal(method,'exact')
    p = min(1, 2 .* binocdf(sgn, n, 0.5));  % p > 1 means center value double-counted
    zval = NaN .* ones(size(p));
else
    % Do a continuity correction, keeping in mind the right direction
    z = (npos-nneg - sign(npos-nneg)) ./ (n .^ .5);
    p = 2 .* normcdf(-abs(z), 0, 1);
    stats.zval = z;
    stats.p = p;
end

% FDR correction
% p = p(stats.isvalid); % already done
pthr = FDR(p, .05);
if isempty(pthr), pthr = Inf; end

stats.pthr_fdr05 = pthr;
sig = p < .05;
fdr_sig = p < pthr;


% Re-insert removed CASES (rows) and fill with zeros
mydirection = zeroinsert(~isvalid, mydirection);
n = zeroinsert(~isvalid, n);
p = zeroinsert(~isvalid, p);
zval = zeroinsert(~isvalid, zval);
sig = zeroinsert(~isvalid, sig);
fdr_sig = zeroinsert(~isvalid, fdr_sig);


stats = struct('method', method, 'isvalid', isvalid, 'n', n, ...
'direction', mydirection, 'p', p, 'zval', zval, 'pthr_fdr05', pthr, ...
'sig', sig, 'fdr_sig', fdr_sig);


end % function


    
