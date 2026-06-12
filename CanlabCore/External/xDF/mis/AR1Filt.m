function ts=AR1Filt(Y,L)
%   ts=AR1Filt(Y,L)
%   Y  should be a vector or a matrix of voxel x time-points
%   L  is number of time-points
%   ts filtered time series
%
% Low/high pass filter for removing AR1.  
% Transfer Function, H(z)=1/(1+AR1*z^-1)
%
%
% SA, Ox, 2017
% srafyouni@gmail.com 
%
%%%REFERENCE
% Arbabshirani et al (2014). Impact of Autocorrelation on Functional 
% Connectivity. NeuroImage. http://doi.org/10.1016/j.neuroimage.2014.07.045
%
if size(Y,2)~=L
    Y=Y';
end

Y=Y-mean(Y,2); 
%only works on >2016 Matlabs, but faster!
% if <2016, use Y=Y-repmat(mean(Y,2),1,L) instead.

%AR1 = autocorr(Y,1);
%AR1 = AR1(2);

AC  = AC_fft(Y,L);
AR1 = AC(:,2);

for v=1:numel(AR1)
    ts(v,:) = filter(1,[1 AR1(v)],Y(v,:));
end