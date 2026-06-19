function [p_unadj,p_adj]=AC_test(xAC,varargin)
% [p_unadj,p_adj]=AC_test(xAC,varargin)
%
% Test the autocorrelation function, following Anderson's a.c.f variance
% estimation. NB! Only T/5 of the lags were used to estimate the variance. 
% p-values can be corrected via 'FDR' or 'Bon' option. 
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
alp = 0.05; 
L   = size(xAC,2);
I   = size(xAC,1);

if round(mean(xAC(:,1)),2)==1
    xAC(:,1)=[];
end

varacf  = (1+2.*sum(xAC(:,1:(L/5)).^2,2))./L; %From Anderson's p8: variance of a.c.f
zs      = xAC./sqrt(varacf);     %z-scores
acpvals = 2.*normcdf(-abs(zs));  %pvals

p_unadj=[ones(I,1)  acpvals];

p_adj=zeros(size(xAC));

if sum(strcmpi(varargin,'FDR')) || nargin==1
    disp('FDR performed.')
    %FDR
    for i=1:I 
        p_adj(i,:) = fdr_bh(acpvals(i,:)); %fdr_bh is an external func: https://uk.mathworks.com/matlabcentral/fileexchange/27418-fdr-bh
    end
elseif sum(strcmpi(varargin,'Bon'))
    disp('Bon performed.')
    p_adj = acpvals;
    p_adj(p_adj<(alp./L)) = 1;
    p_adj(p_adj~=1)       = 0;
end

p_adj=[ones(I,1)  p_adj];