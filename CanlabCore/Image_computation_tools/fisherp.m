function [z,p,sig,pt] = fisherp(p,varargin)
% :Usage:
% ::
%
%     function [z,p,sig,pt] = fisherp(p,[alph])
% 
% :Inputs:
%
%   **p:**
%        values in 4-D array
%
%        1st 3 dims are within images, dim4 = image
%
% :Optional: alpha value for thresholding
% 
% :Outputs:
%
%   **z:**
%        Fisher's combined test statistic, compare to normal
%
%   **p:**
%        p-values for combined test
%
%   **sig:**
%        signficance 1 / 0 binary mask, p < .05 (or alph) FDR-corr
%
%   **pt:**
%        p-value threshold for FDR corrected significance at alph
%
% :Described in:
% Lazar, N. A., Luna, B., Sweeney, J. A., & Eddy, W. F. (2002). 
% Combining brains: a survey of methods for statistical pooling 
% of information. Neuroimage, 16(2), 538-550.
%
% Stouffer, S. A., Suchman, E. A., DeVinney, L. C., Star, S. A., and
% Williams, R. M. 1949. The American Soldier: Vol. I. Adjustment
% During Army Life. Princeton University Press, Princeton.
% 
% Threshold is determined with False Discovery Rate (Benjamini & Hochberg, 1995)
%
% ..
%    tor wager
% ..


if length(varargin) > 0, alph = varargin{1};,else,alph = 0.05;,end

lastdim = length(size(p));
k = size(p,lastdim);

z = -2*sum(log(p),lastdim);
p = 1 - chi2cdf(z(:),2*k);  % distributed as chi-square with 2k df
p = reshape(p,size(z));

% eliminate voxels outside of brain, all 0, 1, or all NaN values
p(all(p == 0,lastdim)) = NaN;
p(all(p == 1,lastdim)) = NaN;

pp = p(:); pp(isnan(pp)) = [];
pt = FDR(pp,alph);
if isempty(pt), pt = -Inf;, end

%z = sum(norminv(1 - p),lastdim) ./ sqrt(k);
%p = normcdf(1-z);

%if alph, sig = p <= alph;,end
sig = p <= pt;

return
