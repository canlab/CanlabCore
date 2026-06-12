function [V,Z,P,BCF] = Bartlett46_fft(Y,T)
% [BCF,BCFA]=Bartlett46_fft(Y,L)
% Super-fats calculation of Bartlett correction factors. 
%
%%%INPUTS:
%   Y: IxT or TxI BOLD voxel-wise/pacellated time series 
%   L: Length of the time series! Just to check the dimensions are all
%   sorted!
%
%%%OUTPUTS:
%   BCF:  In the classic 1935 Bartlett's Correction Factor. \hat{N}=N/BCF
%
%%%DEPENDENCIES:
%   AC_fft.m 
% 
%%%REFERENCES:
%   Bretherton et al, Journal of Climate, 1999, p2004
%   Robert Haining, Geographical Analysis, 1991
%   Richardson & Clifford, 1991, p300
%
%  Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols. 
% "Effective Degrees of Freedom of the Pearson's Correlation Coefficient 
%  under Serial Correlation." bioRxiv (2018): 453795.
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________
if size(Y,1)~=T
    Y=Y'; %IxT
end

xAC = AC_fft(Y,T);
xAC = xAC(:,2); % AC lag-1

xAC = xAC*xAC';

BCF = (1+xAC)./(1-xAC);
V   = BCF./T;

Z   = atanh(corr(Y)).*sqrt(T./BCF-3);
P   = 2 * normcdf(-abs(Z));
%BCFA=[]; %BCFA = (1+xAC).*((1-corr(Y').^2).^2)./(1-xAC);

if any(any(BCF>T))
    BCF (BCF>T) = T;
end