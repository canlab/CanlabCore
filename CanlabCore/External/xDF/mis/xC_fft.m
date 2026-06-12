function [xC,lidx]=xC_fft(Y,T,varargin)
%[xAC]=xC_fft(Y,T,varargin)
%   Super fast full-lag cross-correlation calculation of multi-dimensional 
%   matrices.
%
%%%%%  INPUTS:
%   Y:      A matrix of size IxT comprised of I time series of T length.
%   L:      Time series length
%
%%%%%  OUTPUTS:
%   xC   :  Is a 3D matrix of IxIxT: 
%
%           1) On the diagolans, you get the *auto*correlation. Although note
%              that the result is a symmetric ACF (negatives and positives)
%              autocorrelation lags.
%
%           2) Off diagonals are *cross*correlations between a pair
%           3) The identity IxIx1 is a correlation matrix (i.e. lag-0 structures).
%
%   lidx :  I a vector of lag indexes. Each 1x1xT structure follows these
%           lag indexes.
%%%%%% NOTES:
%   If you need the Pearson's correlation (lag0 xcorr), set lag to 0! Also,
%   diagonal of layer n is the AC of lag n! 
%
%   This function only works for Matlab 2016< and there is no way around
%   it!
%
%   For only a pair time series, this is slower than crosscorr. Only use 
%   this function if you have a lager number of time series 
%
%%%REFERENCES:
%  Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols. 
% "Effective Degrees of Freedom of the Pearson's Correlation Coefficient 
%  under Serial Correlation." bioRxiv (2018): 453795.
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2018
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

if size(Y,2)~=T
    Y=Y'; %IxT
end
I = size(Y,1);

if sum(strcmpi(varargin,'lag'))
    mxL = varargin{find(strcmpi(varargin,'lag'))+1};
    if ~mxL; mxL=1; end; %the user wants the Pearson's Correlation!
else
    mxL = T;
end

Y = Y-mean(Y,2);
%only works on >2016 Matlabs, but faster!
% if <2016, use Y=Y-repmat(mean(Y,2),1,L) instead.

nfft  = 2^nextpow2(2*T-1); %zero-pad me
yfft  = fft(Y,nfft,2);

mxLcc = (mxL-1)*2+1;
xC    = zeros(I,I,mxLcc);
[xx,yy]=find(triu(ones(I),1));

for i=1:numel(xx)
        
    xC0         = ifft(yfft(xx(i),:).*conj(yfft(yy(i),:)));
    xC0         = fliplr([xC0(end-mxL+2:end),xC0(1:mxL)]);
        
    Norm        = sqrt(sum(abs(Y(xx(i),:)).^2)*sum(abs(Y(yy(i),:)).^2));
        
    xC(xx(i),yy(i),:)   = xC0./Norm;
    clear xC0
end

xC = xC + permute(xC,[2 1 3]);

lidx = [-(mxL-1) : (mxL-1)];

