function [xAC,CI,ACOV]=AC_fft(Y,T,varargin)
%[xAC]=AC_fft(Y,T,varargin)
% Super fast full-lag AC calculation of multi-dimention matrices. The
% function exploits fft to estimate the autocorrelations. 
%
%%%%INPUTS
%   Y:      A matrix of size IxT comprised of I time series of T length.
%   T:      Time series length
%   
%   To get a double sided AC, add 'two-sided' as an input.
%
%   NB! The function includes the 0lag zero (i.e. 1) to the output. 
%%%%OUTPUTS
%   xAC:    IxT-1 matrix of full-lag autocorrelations. If 'two-sided'
%           switched, then the matrix is Ix2(T-1).
%   CI :    95% Confidence Intervals of AFC.
%
%%%REFERENCES:
%  Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols. 
% "Effective Degrees of Freedom of the Pearson's Correlation Coefficient 
%  under Serial Correlation." bioRxiv (2018): 453795.
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

if size(Y,2)~=T
    Y=Y';
end

if size(Y,2)~=T
    error('Use IxT or TxI input form.')
end

Y = Y-mean(Y,2); 
%only works on >2016 Matlabs, but faster!
% if <2016, use Y=Y-repmat(mean(Y,2),1,L) instead.

YNorm   = sum(abs(Y).^2,2);
nfft    = 2.^nextpow2(2*T-1); %zero-pad the hell out!
yfft    = fft(Y,nfft,2); %be careful with the dimensions

clear Y %yeah! that bad!

ACOV = real(ifft(yfft.*conj(yfft),[],2));

ACOV = ACOV(:,1:T);

xAC  = ACOV./YNorm; %normalise the COVs

if sum(strcmpi(varargin,'two-sided')) %two sided is just mirrored, AC func is symmetric
   xAC  = [xAC(:,end-T+2:end) ; xAC(:,1:T)];
else
    xAC  = xAC(:,1:T);
end

bnd=(sqrt(2)*erfinv(0.95))./sqrt(T); %assumes normality for AC
CI=[-bnd bnd];

