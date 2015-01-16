function [z,p,sig] = stouffer(p,varargin)
% function [z,p,sig] = stouffer(p,[alph])
% 
% inputs: p values in 4-D array
%         1st 3 dims are within images, dim4 = image
%         optional: alpha value for thresholding
% 
% outputs:
%   z   stouffer's combined test statistic, compare to normal
%   p   p-values for combined test
%   sig signficance 1 / 0 binary mask, if alpha is specified
%
% tor wager
%
% Described in:
% Lazar, N. A., Luna, B., Sweeney, J. A., & Eddy, W. F. (2002). 
% Combining brains: a survey of methods for statistical pooling 
% of information. Neuroimage, 16(2), 538-550.
%
% Stouffer, S. A., Suchman, E. A., DeVinney, L. C., Star, S. A., and
% Williams, R. M. 1949. The American Soldier: Vol. I. Adjustment
% During Army Life. Princeton University Press, Princeton.
% 

if length(varargin) > 0, alph = varargin{1};,else,alph = 0;,end

lastdim = length(size(p));
k = size(p,lastdim);

z = sum(norminv(1 - p),lastdim) ./ sqrt(k);
p = normcdf(1-z);

if alph, sig = p <= alph;,end

return
